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
    #[inline(always)]
    fn get_substitution_score(&self, _pos: (usize, usize), s: u8, r: u8) -> Score {
        if s == r { self.match_score } else { self.mismatch_penalty }
    }
    
    #[inline(always)]
    fn get_substitution_score_v(&self, _pos: (usize, usize), subjects: i16x8, reference: u8) -> i16x8 {
        let v_ref = i16x8::from(reference as i16);
        let v_is_match = subjects.cmp_eq(v_ref);
        v_is_match.blend(i16x8::from(self.match_score), i16x8::from(self.mismatch_penalty))
    }

    #[inline(always)]
    fn get_subject_gap_opening_penalty(&self, _pos: usize) -> Score {
        self.subject_gap_penalty
    }
    #[inline(always)]
    fn get_reference_gap_opening_penalty(&self, _pos: usize) -> Score {
        self.reference_gap_penalty
    }
    #[inline(always)]
    fn get_max_reference_size(&self) -> usize {
        self.max_reference_size
    }
    #[inline(always)]
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
}

impl<C: AlignmentConfig> Aligner<C> for GlobalNtAligner<C> {
    fn align(&mut self, subject: &[u8]) -> Result<Alignment, AlignmentError> {
        let results = self.align_batch(&[subject]);
        results.into_iter().next().expect("align_batch must return exactly one result for a single subject input")
    }

    /// Aligns a batch of subject sequences against the aligner's fixed reference sequence.
    ///
    /// This implementation uses **inter-sequence SIMD vectorization** via the `wide` crate,
    /// allowing 8 independent alignments to be processed simultaneously in a single CPU instruction stream.
    /// It also employs **row recycling** to reduce the score matrix memory footprint from $O(N \times M)$
    /// to $O(M)$, ensuring the "active" scores stay within the L1/L2 caches.
    ///
    /// The algorithm follows the **Needleman-Wunsch** global alignment pattern.
    ///
    /// # Arguments
    /// * `subjects` - A slice of byte slices, each representing a nucleotide sequence to be aligned.
    ///
    /// # Returns
    /// A `Vec` of results, where each item is either an `Alignment` or an `AlignmentError`.
    fn align_batch(&mut self, subjects: &[&[u8]]) -> Vec<Result<Alignment, AlignmentError>> {
        let n_subjects = subjects.len();
        if n_subjects == 0 { return Vec::new(); }

        let mut all_results = Vec::with_capacity(n_subjects);
        let ref_len = self.reference.len();
        let ncols = ref_len + 1;

        // Process subjects in chunks of 8 to match i16x8 SIMD width
        for chunk_idx in (0..n_subjects).step_by(8) {
            let chunk_end = (chunk_idx + 8).min(n_subjects);
            let chunk_subjects = &subjects[chunk_idx..chunk_end];
            let actual_batch_size = chunk_subjects.len();

            let max_sub_len = chunk_subjects.iter().map(|s| s.len()).max().unwrap_or(0);
            let nrows = max_sub_len + 1;

            // Re-use SIMD buffers to avoid expensive heap allocations per batch.
            // We use two rows of i16x8 scores (current and previous) to keep the active
            // working set within the L1/L2 caches (Row Recycling).
            if self.scores[0].len() < ncols {
                self.scores[0].resize(ncols, i16x8::ZERO);
                self.scores[1].resize(ncols, i16x8::ZERO);
            }
            // The ops buffer stores the directions for all 8 alignments in the batch.
            if self.ops.len() < nrows * ncols {
                self.ops.resize(nrows * ncols, i16x8::ZERO);
            }

            // Fill Top Row: Broadcast the precalculated scalar row-0 scores and ops
            // into all 8 lanes of the SIMD buffers.
            for col in 0..ncols {
                self.scores[0][col] = i16x8::from(self.top_row_scores[col]);
                self.ops[col] = i16x8::from(self.top_row_ops[col] as i16);
            }

            let mut final_scores = [0i16; 8];
            
            // Vectorized Fill Phase: Iterate through rows (subject bases).
            // Each cell in the matrix represents the optimal score for the alignment
            // of the first i bases of a subject and the first j bases of the reference.
            for row in 1..nrows {
                let curr = row % 2;
                let prev = (row - 1) % 2;

                // Transpose: Gather the i-th base of all 8 subjects into one SIMD register.
                // This enables inter-sequence parallelism using the Needleman-Wunsch recurrence.
                let mut sub_bases = [0i16; 8];
                for (i, sub) in chunk_subjects.iter().enumerate() {
                    if row <= sub.len() {
                        sub_bases[i] = sub[row - 1] as i16;
                    }
                }
                let v_sub_bases = i16x8::from(sub_bases);

                // Initialize column 0 for this row (Insert state).
                let ref_gap_penalty = self.config.get_reference_gap_opening_penalty(row - 1);
                let v_ref_gap = i16x8::from(ref_gap_penalty);
                let v_acc_ref = self.scores[prev][0] + v_ref_gap;
                self.scores[curr][0] = v_acc_ref;
                self.ops[row * ncols] = i16x8::from(Op::INSERT as i16);

                // Sweep across the reference bases.
                for col in 1..ncols {
                    let r = self.reference[col - 1];

                    // Use the vectorized interface method to support position-specific scores efficiently.
                    let v_sub_score = self.config.get_substitution_score_v((row, col), v_sub_bases, r);

                    // Recurrence: NW(i, j) = max(diag + substitution, up + gap, left + gap)
                    let v_score_match = self.scores[prev][col - 1] + v_sub_score;
                    
                    // Retrieve gap penalties via the interface
                    let v_ref_gap = i16x8::from(self.config.get_reference_gap_opening_penalty(row - 1));
                    let v_sub_gap = i16x8::from(self.config.get_subject_gap_opening_penalty(col - 1));

                    let v_score_insert = self.scores[prev][col] + v_ref_gap;
                    let v_score_delete = self.scores[curr][col - 1] + v_sub_gap;

                    let v_max_score = v_score_match.max(v_score_insert.max(v_score_delete));
                    
                    // Traceback Encoding: Store the winning operation in each lane.
                    let mask_match = v_max_score.cmp_eq(v_score_match);
                    let mask_insert = v_max_score.cmp_eq(v_score_insert) & !mask_match;

                    let v_ops = mask_match.blend(i16x8::from(Op::MATCH as i16), 
                                    mask_insert.blend(i16x8::from(Op::INSERT as i16), 
                                                    i16x8::from(Op::DELETE as i16)));

                    self.scores[curr][col] = v_max_score;
                    self.ops[row * ncols + col] = v_ops;
                }

                // Sequence completion check: If a subject is finished at this row,
                // capture its final alignment score from the last column.
                for (i, sub) in chunk_subjects.iter().enumerate() {
                    if row == sub.len() {
                        let row_scores: [i16; 8] = self.scores[curr][ref_len].into();
                        final_scores[i] = row_scores[i];
                    }
                }
            }

            // Handle edge case: empty subject sequences.
            for (i, sub) in chunk_subjects.iter().enumerate() {
                if sub.len() == 0 {
                    let row0_scores: [i16; 8] = self.scores[0][ref_len].into();
                    final_scores[i] = row0_scores[i];
                }
            }

            // Traceback Phase: Reconstruct the optimal alignment path for each sequence.
            // We read directly from the SIMD 'ops' buffer, extracting one lane at a time.
            for i in 0..actual_batch_size {
                let sub = chunk_subjects[i];
                let end_idx = (sub.len(), ref_len);

                let mut builder = AlignmentBuilder::new(sub, &self.reference);
                let mut cursor = end_idx;
                while cursor != (0, 0) {
                    let l_idx = cursor.0 * ncols + cursor.1;
                    // Extract the Op for the i-th sequence in the batch.
                    let ops_simd: [i16; 8] = self.ops[l_idx].into();
                    let op: Op = unsafe { std::mem::transmute(ops_simd[i] as u8) };
                    builder.take(op, cursor);
                    cursor = matrix::move_back_op(op, cursor);
                }
                builder.take(Op::START, cursor);
                all_results.push(Ok(builder.build(final_scores[i])));
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
