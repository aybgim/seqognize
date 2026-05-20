use crate::config::AlignmentConfig;
use crate::aligner::{Aligner};
use crate::alignment::{Alignment, AlignmentBuilder};
use crate::matrix::{Matrix, Idx, AlignmentError};
use crate::{matrix};
use crate::element::{Score, Element, Op};
use wide::*;

pub struct NtAlignmentConfig {
    pub match_score: Score,
    pub mismatch_penalty: Score,
    pub subject_gap_penalty: Score,
    pub reference_gap_penalty: Score,
    pub max_reference_size: usize,
    pub max_subject_size: usize,
}

impl NtAlignmentConfig {
    pub fn new(match_score: Score, mismatch_penalty: Score, subject_gap_penalty: Score, reference_gap_penalty: Score) -> Self {
        let p_max = match_score.abs()
            .max(mismatch_penalty.abs())
            .max(subject_gap_penalty.abs())
            .max(reference_gap_penalty.abs());
        
        let limit = if p_max > 0 { (16383 / p_max) as usize } else { usize::MAX };
        NtAlignmentConfig {
            match_score,
            mismatch_penalty,
            subject_gap_penalty,
            reference_gap_penalty,
            max_reference_size: limit,
            max_subject_size: limit,
        }
    }
}

impl AlignmentConfig for NtAlignmentConfig {
    fn get_substitution_score(&self, _pos: (usize, usize), s: u8, r: u8) -> Score {
        if s == r { self.match_score } else { self.mismatch_penalty }
    }
    fn get_subject_gap_opening_penalty(&self, _pos: usize) -> Score {
        self.subject_gap_penalty
    }
    fn get_reference_gap_opening_penalty(&self, _pos: usize) -> Score {
        self.reference_gap_penalty
    }
    fn get_max_reference_size(&self) -> usize {
        self.max_reference_size
    }
    fn get_max_subject_size(&self) -> usize {
        self.max_subject_size
    }
}

pub struct GlobalNtAligner {
    pub config: NtAlignmentConfig,
    pub reference: Vec<u8>,
    pub top_row_scores: Vec<Score>,
    pub top_row_ops: Vec<Op>,
    pub scores: [Vec<i16x8>; 2],
    pub ops: Vec<i16x8>,
}

impl GlobalNtAligner {
    pub fn new(config: NtAlignmentConfig, reference: Vec<u8>) -> Self {
        let ncols = reference.len() + 1;
        let mut top_row_scores = vec![0; ncols];
        let mut top_row_ops = vec![Op::START; ncols];
        
        // Only pre-calculate if within safe limits to avoid i16 overflow
        if reference.len() <= config.get_max_reference_size() {
            let mut acc = 0;
            for col in 1..ncols {
                acc += config.get_subject_gap_opening_penalty(col - 1);
                top_row_scores[col] = acc;
                top_row_ops[col] = Op::DELETE;
            }
        }

        GlobalNtAligner {
            config,
            reference,
            top_row_scores,
            top_row_ops,
            scores: [Vec::new(), Vec::new()],
            ops: Vec::new(),
        }
    }
}

impl Aligner<NtAlignmentConfig> for GlobalNtAligner {
    fn align(&mut self, subject: &[u8]) -> Result<Alignment, AlignmentError> {
        let results = self.align_batch(&[subject]);
        results.into_iter().next().unwrap()
    }

    fn align_batch(&mut self, subjects: &[&[u8]]) -> Vec<Result<Alignment, AlignmentError>> {
        let n_subjects = subjects.len();
        if n_subjects == 0 { return Vec::new(); }

        let mut all_results = Vec::with_capacity(n_subjects);
        let ref_len = self.reference.len();
        let ncols = ref_len + 1;

        for chunk_idx in (0..n_subjects).step_by(8) {
            let chunk_end = (chunk_idx + 8).min(n_subjects);
            let chunk_subjects = &subjects[chunk_idx..chunk_end];
            let actual_batch_size = chunk_subjects.len();

            let max_sub_len = chunk_subjects.iter().map(|s| s.len()).max().unwrap_or(0);
            let nrows = max_sub_len + 1;

            // Prepare SIMD buffers
            if self.scores[0].len() < ncols {
                self.scores[0].resize(ncols, i16x8::ZERO);
                self.scores[1].resize(ncols, i16x8::ZERO);
            }
            if self.ops.len() < nrows * ncols {
                self.ops.resize(nrows * ncols, i16x8::ZERO);
            }

            // Fill Top Row (broadcast precalculated row 0)
            for col in 0..ncols {
                self.scores[0][col] = i16x8::from(self.top_row_scores[col]);
                self.ops[col] = i16x8::from(self.top_row_ops[col] as i16);
            }

            let match_score_v = i16x8::from(self.config.match_score);
            let mismatch_penalty_v = i16x8::from(self.config.mismatch_penalty);
            let ref_gap_v = i16x8::from(self.config.reference_gap_penalty);
            let sub_gap_v = i16x8::from(self.config.subject_gap_penalty);

            // Active masks for varying lengths
            let mut active_masks = [i16x8::ZERO; 1]; // We don't strictly need a mask if we just align dummy data, but good for safety.
            
            for row in 1..nrows {
                let curr = row % 2;
                let prev = (row - 1) % 2;

                // Load subject bases for this row across all 8 subjects
                let mut sub_bases = [0i16; 8];
                for (i, sub) in chunk_subjects.iter().enumerate() {
                    if row <= sub.len() {
                        sub_bases[i] = sub[row - 1] as i16;
                    }
                }
                let v_sub_bases = i16x8::from(sub_bases);

                // Initialize column 0 for this row
                let v_acc_ref = self.scores[prev][0] + ref_gap_v;
                self.scores[curr][0] = v_acc_ref;
                self.ops[row * ncols] = i16x8::from(Op::INSERT as i16);

                for col in 1..ncols {
                    let v_ref_base = i16x8::from(self.reference[col - 1] as i16);
                    let v_is_match = v_sub_bases.cmp_eq(v_ref_base);
                    let v_sub_score = v_is_match.blend(match_score_v, mismatch_penalty_v);

                    let v_score_match = self.scores[prev][col - 1] + v_sub_score;
                    let v_score_insert = self.scores[prev][col] + ref_gap_v;
                    let v_score_delete = self.scores[curr][col - 1] + sub_gap_v;

                    let v_max_score = v_score_match.max(v_score_insert.max(v_score_delete));
                    
                    let mask_match = v_max_score.cmp_eq(v_score_match);
                    let mask_insert = v_max_score.cmp_eq(v_score_insert) & !mask_match;

                    let v_ops = mask_match.blend(i16x8::from(Op::MATCH as i16), 
                                    mask_insert.blend(i16x8::from(Op::INSERT as i16), 
                                                    i16x8::from(Op::DELETE as i16)));

                    self.scores[curr][col] = v_max_score;
                    self.ops[row * ncols + col] = v_ops;
                }
            }

            // Extract results and perform traceback
            for i in 0..actual_batch_size {
                let sub = chunk_subjects[i];
                let end_idx = (sub.len(), ref_len);
                
                // Build a dummy matrix for traceback (transposing from i16x8 back to standard layout)
                // This is suboptimal but allows reusing existing traceback logic.
                // For antibodies, traceback is tiny compared to fill.
                let mut mtx = Matrix::of(sub.len() + 1, ref_len + 1);
                for r in 0..=sub.len() {
                    for c in 0..=ref_len {
                        let l_idx = r * ncols + c;
                        let ops_simd: [i16; 8] = self.ops[l_idx].into();
                        mtx.set_op(r, c, unsafe { std::mem::transmute(ops_simd[i] as u8) });
                    }
                }
                let final_scores_simd: [i16; 8] = self.scores[sub.len() % 2][ref_len].into();
                let final_score = final_scores_simd[i];

                let mut builder = AlignmentBuilder::new(sub, &self.reference);
                let mut cursor = end_idx;
                while cursor != (0, 0) {
                    let op = mtx.get_op(cursor.0, cursor.1);
                    builder.take(op, cursor);
                    cursor = matrix::move_back_op(op, cursor);
                }
                builder.take(Op::START, cursor);
                all_results.push(Ok(builder.build(final_score)));
            }
        }
        all_results
    }

    fn reference(&self) -> &[u8] {
        &self.reference
    }

    fn check_sizes(&self, subject_len: usize, reference_len: usize) -> Result<(), AlignmentError> {
        if subject_len > self.config.get_max_subject_size() || reference_len > self.config.get_max_reference_size() {
            return Err(AlignmentError::SequenceTooLong);
        }
        Ok(())
    }

    fn fill_top_row(&self, _mtx: &mut Matrix) {}
    fn fill_left_column(&self, _mtx: &mut Matrix) {}
    fn fill(&self, _mtx: &mut Matrix, _subject: &[u8], _reference: &[u8]) {}
    fn end_idx(&self, _mtx: &Matrix) -> Idx { (0, 0) }
    fn trace_back(&self, _mtx: &Matrix, _idx: Idx, _s: &[u8], _r: &[u8]) -> Alignment { 
        Alignment::from("", "", 0) 
    }
}

pub fn insertion(score: Score) -> Element {
    Element { op: Op::INSERT, score }
}

pub fn deletion(score: Score) -> Element {
    Element { op: Op::DELETE, score }
}

pub fn substitution(score: Score) -> Element {
    Element { op: Op::MATCH, score }
}

#[cfg(test)]
mod tests {
    use crate::nt_aligner::{GlobalNtAligner, NtAlignmentConfig, deletion, substitution};
    use crate::aligner::Aligner;
    use crate::matrix::AlignmentError;
    use crate::alignment::Alignment;
    use crate::element::{Score, Element, Op};

    fn aligner(reference: &[u8]) -> GlobalNtAligner {
        GlobalNtAligner::new(
            NtAlignmentConfig::new(1, -1, -1, -1),
            reference.to_vec()
        )
    }

    #[test]
    fn test_match() {
        assert_eq!(
            aligner(b"AGCT").align(b"AGCT").unwrap(),
            Alignment::from("AGCT", "AGCT", 4)
        )
    }

    #[test]
    fn test_mismatch() {
        assert_eq!(
            aligner(b"AGCT").align(b"AGAT").unwrap(),
            Alignment::from("AGAT", "AGCT", 2)
        )
    }

    #[test]
    fn test_double_insertion() {
        assert_eq!(
            aligner(b"AT").align(b"AGCT").unwrap(),
            Alignment::from("AGCT", "A__T", 0)
        )
    }

    #[test]
    fn test_batch_alignment() {
        let mut al = aligner(b"AGCT");
        let subjects = vec![b"AGCT".as_slice(), b"AGAT".as_slice(), b"AG_T".as_slice()];
        let results = al.align_batch(&subjects);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].as_ref().unwrap().score, 4);
        assert_eq!(results[1].as_ref().unwrap().score, 2);
    }
}
