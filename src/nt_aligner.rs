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
    
    #[inline(always)]
    fn get_substitution_score_v(&self, _pos: (usize, usize), subjects: i16x8, reference: u8) -> i16x8 {
        let v_ref = i16x8::from(reference as i16);
        let v_is_match = subjects.cmp_eq(v_ref);
        v_is_match.blend(i16x8::from(self.match_score), i16x8::from(self.mismatch_penalty))
    }

    fn get_subject_gap_opening_penalty(&self, _pos: usize) -> Score {
        self.subject_gap_penalty
    }
    #[inline(always)]
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

pub struct GlobalNtAligner<C: AlignmentConfig> {
    pub config: C,
    pub reference: Vec<u8>,
    pub top_row_scores: Vec<Score>,
    pub top_row_ops: Vec<Op>,
    pub scores: [Vec<i16x8>; 2],
    pub ops: Vec<i16x8>,
}

impl<C: AlignmentConfig> GlobalNtAligner<C> {
    pub fn new(config: C, reference: Vec<u8>) -> Self {
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

    #[inline]
    fn prepare_batch_buffers(&mut self, nrows: usize, ncols: usize) {
        if self.scores[0].len() < ncols {
            self.scores[0].resize(ncols, i16x8::ZERO);
            self.scores[1].resize(ncols, i16x8::ZERO);
        }
        if self.ops.len() < nrows * ncols {
            self.ops.resize(nrows * ncols, i16x8::ZERO);
        }
    }

    #[inline]
    fn initialize_fill(&mut self, ncols: usize) {
        for col in 0..ncols {
            self.scores[0][col] = i16x8::from(self.top_row_scores[col]);
            self.ops[col] = i16x8::from(self.top_row_ops[col] as i16);
        }
    }

    #[inline]
    fn compute_fill(&mut self, chunk_subjects: &[&[u8]], nrows: usize, ncols: usize) -> [i16; 8] {
        let ref_len = self.reference.len();
        let mut final_scores = [0i16; 8];

        for row in 1..nrows {
            let curr = row % 2;
            let prev = (row - 1) % 2;

            let mut sub_bases = [0i16; 8];
            for (i, sub) in chunk_subjects.iter().enumerate() {
                if row <= sub.len() {
                    sub_bases[i] = sub[row - 1] as i16;
                }
            }
            let v_sub_bases = i16x8::from(sub_bases);

            let ref_gap_penalty = self.config.get_reference_gap_opening_penalty(row - 1);
            let v_ref_gap = i16x8::from(ref_gap_penalty);
            self.scores[curr][0] = self.scores[prev][0] + v_ref_gap;
            self.ops[row * ncols] = i16x8::from(Op::INSERT as i16);

            for col in 1..ncols {
                let r = self.reference[col - 1];
                let v_sub_score = self.config.get_substitution_score_v((row, col), v_sub_bases, r);

                let v_score_match = self.scores[prev][col - 1] + v_sub_score;
                let v_ref_gap = i16x8::from(self.config.get_reference_gap_opening_penalty(row - 1));
                let v_sub_gap = i16x8::from(self.config.get_subject_gap_opening_penalty(col - 1));

                let v_score_insert = self.scores[prev][col] + v_ref_gap;
                let v_score_delete = self.scores[curr][col - 1] + v_sub_gap;

                let v_max_score = v_score_match.max(v_score_insert.max(v_score_delete));
                
                let mask_match = v_max_score.cmp_eq(v_score_match);
                let mask_insert = v_max_score.cmp_eq(v_score_insert) & !mask_match;

                let v_ops = mask_match.blend(i16x8::from(Op::MATCH as i16), 
                                mask_insert.blend(i16x8::from(Op::INSERT as i16), 
                                                i16x8::from(Op::DELETE as i16)));

                self.scores[curr][col] = v_max_score;
                self.ops[row * ncols + col] = v_ops;
            }

            for (i, sub) in chunk_subjects.iter().enumerate() {
                if row == sub.len() {
                    let row_scores: [i16; 8] = self.scores[curr][ref_len].into();
                    final_scores[i] = row_scores[i];
                }
            }
        }

        for (i, sub) in chunk_subjects.iter().enumerate() {
            if sub.len() == 0 {
                let row0_scores: [i16; 8] = self.scores[0][ref_len].into();
                final_scores[i] = row0_scores[i];
            }
        }
        final_scores
    }

    #[inline]
    fn perform_tracebacks(&self, chunk_subjects: &[&[u8]], final_scores: [i16; 8], ncols: usize, all_results: &mut Vec<Result<Alignment, AlignmentError>>) {
        let ref_len = self.reference.len();
        for i in 0..chunk_subjects.len() {
            let sub = chunk_subjects[i];
            let mut builder = AlignmentBuilder::new(sub, &self.reference);
            let mut cursor = (sub.len(), ref_len);
            while cursor != (0, 0) {
                let l_idx = cursor.0 * ncols + cursor.1;
                let ops_simd: [i16; 8] = self.ops[l_idx].into();
                let op: Op = unsafe { std::mem::transmute(ops_simd[i] as u8) };
                builder.take(op, cursor);
                cursor = matrix::move_back_op(op, cursor);
            }
            builder.take(Op::START, cursor);
            all_results.push(Ok(builder.build(final_scores[i])));
        }
    }
}

impl<C: AlignmentConfig> Aligner<C> for GlobalNtAligner<C> {
    fn align(&mut self, subject: &[u8]) -> Result<Alignment, AlignmentError> {
        let results = self.align_batch(&[subject]);
        results.into_iter().next().expect("align_batch must return exactly one result for a single subject input")
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

            let max_sub_len = chunk_subjects.iter().map(|s| s.len()).max().unwrap_or(0);
            let nrows = max_sub_len + 1;

            self.prepare_batch_buffers(nrows, ncols);
            self.initialize_fill(ncols);
            let final_scores = self.compute_fill(chunk_subjects, nrows, ncols);
            self.perform_tracebacks(chunk_subjects, final_scores, ncols, &mut all_results);
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
    use crate::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
    use crate::aligner::Aligner;
    use crate::alignment::Alignment;

    fn aligner(reference: &[u8]) -> GlobalNtAligner<NtAlignmentConfig> {
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
        assert_eq!(results[0].as_ref().expect("Alignment 0 failed").score, 4);
        assert_eq!(results[1].as_ref().expect("Alignment 1 failed").score, 2);
    }
}
