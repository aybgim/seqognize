use crate::aligner::{Aligner, AlignmentError};
use crate::alignment::{Alignment, AlignmentBuilder, Idx};
use crate::config::{AlignmentConfig, Score, SimdScore, LANES};
use crate::alignment::Op;
use crate::aligner::AlignmentError::SequenceTooLong;

use std::convert::TryFrom;
use wide::CmpEq;

/// Configuration for nucleotide alignment scoring.
pub struct NtAlignmentConfig {
    /// Score awarded for matching nucleotide bases.
    pub match_score: Score,
    /// Penalty for mismatching nucleotide bases.
    pub mismatch_penalty: Score,
    /// Penalty for opening a gap in the subject sequence.
    pub subject_gap_penalty: Score,
    /// Penalty for opening a gap in the reference sequence.
    pub reference_gap_penalty: Score,
    /// Maximum allowed reference sequence size.
    pub max_reference_size: usize,
    /// Maximum allowed subject sequence size.
    pub max_subject_size: usize,
}

impl NtAlignmentConfig {
    /// Creates a new `NtAlignmentConfig` with the specified scoring parameters.
    ///
    /// # Arguments
    /// * `match_score` - The reward for matching nucleotide bases.
    /// * `mismatch_penalty` - The penalty for mismatching nucleotide bases (usually negative).
    /// * `subject_gap_penalty` - The penalty for opening a gap in the subject sequence (usually negative).
    /// * `reference_gap_penalty` - The penalty for opening a gap in the reference sequence (usually negative).
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
    /// Returns the substitution score for a pair of nucleotide bases.
    ///
    /// # Arguments
    /// * `_pos` - The (row, col) position in the matrix (unused in this implementation).
    /// * `s` - The subject nucleotide base.
    /// * `r` - The reference nucleotide base.
    fn get_substitution_score(&self, _pos: (usize, usize), s: u8, r: u8) -> Score {
        if s == r { self.match_score } else { self.mismatch_penalty }
    }
    
    /// Returns a SIMD vector of substitution scores for a batch of subject bases against a single reference base.
    ///
    /// # Arguments
    /// * `_pos` - The (row, col) position in the matrix (unused in this implementation).
    /// * `subjects` - A SIMD vector containing subject nucleotide bases.
    /// * `reference` - The reference nucleotide base.
    #[inline(always)]
    fn get_substitution_score_v(&self, _pos: (usize, usize), subjects: SimdScore, reference: u8) -> SimdScore {
        let v_ref = SimdScore::from(reference as i16);
        let v_is_match = subjects.cmp_eq(v_ref);
        v_is_match.blend(SimdScore::from(self.match_score), SimdScore::from(self.mismatch_penalty))
    }

    /// Returns the gap opening penalty for the subject sequence.
    ///
    /// # Arguments
    /// * `_pos` - The position in the sequence (unused in this implementation).
    fn get_subject_gap_opening_penalty(&self, _pos: usize) -> Score {
        self.subject_gap_penalty
    }

    /// Returns the gap opening penalty for the reference sequence.
    ///
    /// # Arguments
    /// * `_pos` - The position in the sequence (unused in this implementation).
    #[inline(always)]
    fn get_reference_gap_opening_penalty(&self, _pos: usize) -> Score {
        self.reference_gap_penalty
    }

    /// Returns the maximum allowed size for the reference sequence.
    fn get_max_reference_size(&self) -> usize {
        self.max_reference_size
    }

    /// Returns the maximum allowed size for the subject sequence.
    fn get_max_subject_size(&self) -> usize {
        self.max_subject_size
    }
}

/// A global nucleotide aligner that uses SIMD-accelerated Needleman-Wunsch.
pub struct GlobalNtAligner<C: AlignmentConfig> {
    /// The alignment configuration (scoring, size limits).
    pub config: C,
    /// The fixed reference sequence to align against.
    pub reference: Vec<u8>,
    /// Precomputed scores for the first row of the dynamic programming matrix.
    pub top_row_scores: Vec<Score>,
    /// Precomputed operations for the first row of the dynamic programming matrix.
    pub top_row_ops: Vec<Op>,
    /// Rolling buffers for the dynamic programming matrix scores (only 2 rows needed).
    pub scores: [Vec<SimdScore>; 2],
    /// Compressed operation matrix for traceback, stored as SIMD vectors.
    pub ops: Vec<SimdScore>,
}

impl<C: AlignmentConfig> GlobalNtAligner<C> {
    /// Creates a new `GlobalNtAligner` with the given configuration and reference sequence.
    ///
    /// # Arguments
    /// * `config` - The alignment configuration (scoring, size limits).
    /// * `reference` - The fixed reference sequence to align against.
    ///
    /// # Returns
    /// `Ok(GlobalNtAligner)` if successful, or `Err(AlignmentError::SequenceTooLong)` if the reference sequence exceeds the maximum allowed size.
    pub fn new(config: C, reference: Vec<u8>) -> Result<Self, AlignmentError> {
        if reference.len() > config.get_max_reference_size() {
            return Err(AlignmentError::SequenceTooLong);
        }

        let ncols = reference.len() + 1;
        let mut top_row_scores = vec![0; ncols];
        let mut top_row_ops = vec![Op::START; ncols];

        let mut acc = 0;
        for col in 1..ncols {
            acc += config.get_subject_gap_opening_penalty(col - 1);
            top_row_scores[col] = acc;
            top_row_ops[col] = Op::DELETE;
        }

        Ok(GlobalNtAligner {
            config,
            reference,
            top_row_scores,
            top_row_ops,
            scores: [
                vec![SimdScore::ZERO; ncols],
                vec![SimdScore::ZERO; ncols],
            ],
            ops: Vec::new(),
        })
    }
}

impl<C: AlignmentConfig> Aligner<C> for GlobalNtAligner<C> {
    /// Aligns a single subject sequence against the reference.
    ///
    /// # Arguments
    /// * `subject` - A slice of nucleotide bases to align against the reference.
    ///
    /// # Returns
    /// `Ok(Alignment)` containing the result, or `Err(AlignmentError::SequenceTooLong)` if the subject exceeds the limit.
    fn align(&mut self, subject: &[u8]) -> Result<Alignment, AlignmentError> {
        let mut results = self.align_batch(&[subject])?;
        Ok(results.pop().expect("align_batch must return exactly one result for a single subject input"))
    }

    /// Aligns a batch of subject sequences against the aligner's fixed reference sequence.
    ///
    /// This method orchestrates the high-performance alignment pipeline by:
    /// 1. Validating that all subjects are within the maximum allowed length.
    /// 2. Grouping subjects into chunks that match the CPU's SIMD width (8 or 16).
    /// 3. Reusing stateful memory buffers to eliminate heap allocation overhead.
    /// 4. Executing a vectorized Needleman-Wunsch fill phase for each chunk.
    /// 5. Performing individual scalar tracebacks to reconstruct the optimal paths.
    ///
    /// # Arguments
    /// * `subjects` - A slice of nucleotide sequences to be aligned against the reference.
    ///
    /// # Returns
    /// `Ok(Vec<Alignment>)` containing results for each input sequence, or `Err(AlignmentError::SequenceTooLong)` if any sequence exceeds the limit.
    fn align_batch(&mut self, subjects: &[&[u8]]) -> Result<Vec<Alignment>, AlignmentError> {
        let n_subjects = subjects.len();
        if n_subjects == 0 { return Ok(Vec::new()); }

        // Early validation: Check if any subject is too long before allocating or computing.
        let max_allowed = self.config.get_max_subject_size();
        for sub in subjects {
            if sub.len() > max_allowed {
                return Err(SequenceTooLong);
            }
        }

        let mut all_results = Vec::with_capacity(n_subjects);
        let ref_len = self.reference.len();
        let ncols = ref_len + 1;

        for chunk_idx in (0..n_subjects).step_by(LANES) {
            let chunk_end = (chunk_idx + LANES).min(n_subjects);
            let chunk_subjects = &subjects[chunk_idx..chunk_end];

            let max_sub_len = chunk_subjects.iter().map(|s| s.len()).max().unwrap_or(0);
            let nrows = max_sub_len + 1;

            self.prepare_batch_buffers(nrows, ncols);
            self.fill_first_row(ncols);
            let final_scores = self.fill_matrix(chunk_subjects, nrows, ncols);
            self.perform_tracebacks(chunk_subjects, final_scores, ncols, &mut all_results);
        }
        Ok(all_results)
    }
}

impl<C: AlignmentConfig> GlobalNtAligner<C> {
    /// Resizes internal buffers to accommodate the required number of rows and columns.
    ///
    /// Reuses existing allocations to minimize heap churn during batch processing.
    ///
    /// # Arguments
    /// * `nrows` - The required number of rows (maximum subject length in the current batch + 1).
    /// * `ncols` - The required number of columns (reference length + 1).
    fn prepare_batch_buffers(&mut self, nrows: usize, ncols: usize) {
        if self.scores[0].len() < ncols {
            self.scores[0].resize(ncols, SimdScore::ZERO);
            self.scores[1].resize(ncols, SimdScore::ZERO);
        }
        
        let unrolled_op_mtx_len = nrows * ncols;
        if self.ops.len() < unrolled_op_mtx_len {
            self.ops.resize(unrolled_op_mtx_len, SimdScore::ZERO);
        }
    }

    /// Initializes the first row of the dynamic programming matrix.
    ///
    /// # Arguments
    /// * `ncols` - The number of columns to initialize.
    fn fill_first_row(&mut self, ncols: usize) {
        for col in 0..ncols {
            self.scores[0][col] = SimdScore::from(self.top_row_scores[col]);
            self.ops[col] = SimdScore::from(self.top_row_ops[col] as i16);
        }
    }

    /// Executes the vectorized Needleman-Wunsch fill phase for a batch of sequences.
    ///
    /// # Arguments
    /// * `chunk_subjects` - The batch of subject sequences to align.
    /// * `nrows` - The maximum number of rows to compute.
    /// * `ncols` - The number of columns in the dynamic programming matrix.
    ///
    /// # Returns
    /// An array containing the final alignment scores for each lane in the SIMD vector.
    fn fill_matrix(&mut self, chunk_subjects: &[&[u8]], nrows: usize, ncols: usize) -> [i16; LANES] {
        let mut final_scores = [0i16; LANES];

        // Row 0 is initialized but not computed in the SIMD loop. We capture 
        // empty subject scores here before the rolling buffer (size 2) 
        // overwrites Row 0 (e.g., Row 2, Row 4, etc. reuse scores[0]).
        self.handle_empty_subjects(chunk_subjects, &mut final_scores);

        for row in 1..nrows {
            let v_sub_bases = self.gather_subject_bases(chunk_subjects, row);
            let v_ref_gap = SimdScore::from(self.config.get_reference_gap_opening_penalty(row - 1));
            self.compute_first_col(row, v_ref_gap, ncols);
            for col in 1..ncols {
                self.compute_cell_simd(row, col, v_sub_bases, v_ref_gap, ncols);
            }
            self.capture_finished_scores(chunk_subjects, row, &mut final_scores);
        }

        final_scores
    }

    /// Gathers nucleotide bases from the current row for all subjects in the batch into a SIMD vector.
    ///
    /// # Arguments
    /// * `chunk_subjects` - The batch of subject sequences.
    /// * `row` - The 1-based index of the current row being processed.
    ///
    /// # Returns
    /// A SIMD vector containing the bases at `row - 1` for each subject, or 0 if the subject is shorter.
    #[inline(always)]
    fn gather_subject_bases(&self, chunk_subjects: &[&[u8]], row: usize) -> SimdScore {
        let mut sub_bases = [0i16; LANES];
        for (i, sub) in chunk_subjects.iter().enumerate() {
            if row <= sub.len() {
                sub_bases[i] = sub[row - 1] as i16;
            }
        }
        SimdScore::from(sub_bases)
    }

    /// Initialize column 0 for this row.
    ///
    /// # Arguments
    /// * `row` - The 1-based index of the current row.
    /// * `v_ref_gap` - A SIMD vector containing the reference gap penalty.
    /// * `ncols` - The total number of columns in the matrix.
    #[inline(always)]
    fn compute_first_col(&mut self, row: usize, v_ref_gap: SimdScore, ncols: usize) {
        let curr = row % 2;
        let prev = (row - 1) % 2;
        self.scores[curr][0] = self.scores[prev][0] + v_ref_gap;
        self.ops[row * ncols] = SimdScore::from(Op::INSERT as i16);
    }

    /// Performs vectorized cell computation for a specific row and column.
    ///
    /// # Arguments
    /// * `row` - The 1-based index of the current row.
    /// * `col` - The 1-based index of the current column.
    /// * `v_sub_bases` - A SIMD vector of nucleotide bases for each subject at this row.
    /// * `v_ref_gap` - A SIMD vector containing the reference gap penalty.
    /// * `ncols` - The total number of columns in the matrix.
    #[inline(always)]
    fn compute_cell_simd(&mut self, row: usize, col: usize, v_sub_bases: SimdScore, v_ref_gap: SimdScore, ncols: usize) {
        let curr = row % 2;
        let prev = (row - 1) % 2;
        let r = self.reference[col - 1];

        // Use the vectorized interface method to support position-specific scores efficiently.
        let v_sub_score = self.config.get_substitution_score_v((row, col), v_sub_bases, r);

        // Recurrence: NW(i, j) = max(diag + substitution, up + gap, left + gap)
        let v_score_match = self.scores[prev][col - 1] + v_sub_score;
        let v_sub_gap = SimdScore::from(self.config.get_subject_gap_opening_penalty(col - 1));

        let v_score_insert = self.scores[prev][col] + v_ref_gap;
        let v_score_delete = self.scores[curr][col - 1] + v_sub_gap;

        let v_max_score = v_score_match.max(v_score_insert).max(v_score_delete);
        
        // Traceback Encoding: Store the winning operation in each lane.
        let mask_match = v_max_score.cmp_eq(v_score_match);
        let mask_insert = v_max_score.cmp_eq(v_score_insert) & !mask_match;

        // Option A: Clean alignment formatting
        let v_ops = mask_match.blend(
            SimdScore::from(Op::MATCH as i16),
            mask_insert.blend(
                SimdScore::from(Op::INSERT as i16),
                SimdScore::from(Op::DELETE as i16),
            ),
        );

        self.scores[curr][col] = v_max_score;
        self.ops[row * ncols + col] = v_ops;
    }

    /// Captures the scores for sequences that have reached their full length at the current row.
    ///
    /// # Arguments
    /// * `chunk_subjects` - The batch of subject sequences.
    /// * `row` - The current row index.
    /// * `final_scores` - A mutable array to store the final scores for each lane.
    #[inline(always)]
    fn capture_finished_scores(&self, chunk_subjects: &[&[u8]], row: usize, final_scores: &mut [i16; LANES]) {
        let curr = row % 2;
        let ref_len = self.reference.len();
        for (i, sub) in chunk_subjects.iter().enumerate() {
            if row == sub.len() {
                let row_scores: [i16; LANES] = self.scores[curr][ref_len].into();
                final_scores[i] = row_scores[i];
            }
        }
    }

    /// Handles sequences with zero length to ensure they get correct gap-only alignment scores.
    ///
    /// # Arguments
    /// * `chunk_subjects` - The batch of subject sequences.
    /// * `final_scores` - A mutable array to store the final scores for each lane.
    #[inline(always)]
    fn handle_empty_subjects(&self, chunk_subjects: &[&[u8]], final_scores: &mut [i16; LANES]) {
        let ref_len = self.reference.len();
        for (i, sub) in chunk_subjects.iter().enumerate() {
            if sub.len() == 0 {
                let row0_scores: [i16; LANES] = self.scores[0][ref_len].into();
                final_scores[i] = row0_scores[i];
            }
        }
    }

    /// Performs scalar traceback for each sequence in the batch to reconstruct the alignment paths.
    ///
    /// # Arguments
    /// * `chunk_subjects` - The batch of subject sequences.
    /// * `final_scores` - The final alignment scores for each lane in the batch.
    /// * `ncols` - The number of columns in the dynamic programming matrix.
    /// * `all_results` - A mutable vector where the reconstructed `Alignment` objects will be stored.
    fn perform_tracebacks(&self, chunk_subjects: &[&[u8]], final_scores: [i16; LANES], ncols: usize, all_results: &mut Vec<Alignment>) {
        let ref_len = self.reference.len();
        for i in 0..chunk_subjects.len() {
            let sub = chunk_subjects[i];
            let mut builder = AlignmentBuilder::new(sub, &self.reference);
            let mut cursor = Idx(sub.len(), ref_len);
            while cursor != Idx(0, 0) {
                let l_idx = cursor.0 * ncols + cursor.1;
                let op = self.to_op(i, l_idx);
                builder.take(op, cursor);
                cursor = cursor.move_back(op);
            }
            builder.take(Op::START, cursor);
            let alignment = builder.build(final_scores[i]);
            all_results.push(alignment);
        }
    }

    /// Converts a SIMD operation vector at the given linear index into a scalar `Op` for a specific lane.
    ///
    /// # Arguments
    /// * `i` - The index of the lane (0 to LANES-1).
    /// * `l_idx` - The linear index into the `ops` buffer.
    ///
    /// # Returns
    /// The `Op` value for the specified lane and position.
    #[inline(always)]
    fn to_op(&self, i: usize, l_idx: usize) -> Op {
        let ops_simd: SimdScore = self.ops[l_idx].into();
        Op::try_from(ops_simd.to_array()[i] as u8).expect("Invalid Op byte!")
    }
}

#[cfg(test)]
mod tests {
    use crate::aligner::{Aligner, AlignmentError};
    use crate::alignment::Alignment;
    use crate::config::Score;
    use crate::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};

    fn aligner(reference: &[u8]) -> GlobalNtAligner<NtAlignmentConfig> {
        GlobalNtAligner::new(
            NtAlignmentConfig::new(1, -1, -1, -1),
            reference.to_vec()
        ).unwrap()
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
    fn test_insertion() {
        assert_eq!(
            aligner(b"AGT").align(b"AGCT").unwrap(),
            Alignment::from("AGCT", "AG_T", 2)
        )
    }

    #[test]
    fn test_deletion() {
        assert_eq!(
            aligner(b"AGCT").align(b"AGT").unwrap(),
            Alignment::from("AG_T", "AGCT", 2)
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
    fn test_double_deletion() {
        assert_eq!(
            aligner(b"AGCT").align(b"AT").unwrap(),
            Alignment::from("A__T", "AGCT", 0)
        )
    }

    #[test]
    fn test_leading_insertion() {
        assert_eq!(
            aligner(b"GCT").align(b"AGCT").unwrap(),
            Alignment::from("AGCT", "_GCT", 2)
        )
    }

    #[test]
    fn test_leading_deletion() {
        assert_eq!(
            aligner(b"AGCT").align(b"GCT").unwrap(),
            Alignment::from("_GCT", "AGCT", 2)
        )
    }

    #[test]
    fn test_trailing_insertion() {
        assert_eq!(
            aligner(b"AGC").align(b"AGCT").unwrap(),
            Alignment::from("AGCT", "AGC_", 2)
        )
    }

    #[test]
    fn test_trailing_deletion() {
        assert_eq!(
            aligner(b"AGCT").align(b"AGC").unwrap(),
            Alignment::from("AGC_", "AGCT", 2)
        )
    }

    #[test]
    fn test_two_insertions() {
        assert_eq!(
            aligner(b"GT").align(b"AGCT").unwrap(),
            Alignment::from("AGCT", "_G_T", 0)
        )
    }

    #[test]
    fn test_two_deletions() {
        assert_eq!(
            aligner(b"AGCT").align(b"AC").unwrap(),
            Alignment::from("A_C_", "AGCT", 0)
        )
    }

    #[test]
    fn test_empty_subject() {
        assert_eq!(
            aligner(b"AGCT").align(b"").unwrap(),
            Alignment::from("____", "AGCT", -4)
        )
    }

    #[test]
    fn test_empty_reference() {
        assert_eq!(
            aligner(b"").align(b"AGCT").unwrap(),
            Alignment::from("AGCT", "____", -4)
        )
    }

    #[test]
    fn test_oversize_subject() {
        let long_seq = vec![b'A'; 40000];
        let result = aligner(b"A").align(&long_seq);
        assert_eq!(result, Err(AlignmentError::SequenceTooLong));
    }

    #[test]
    fn test_oversize_reference() {
        let long_seq = vec![b'A'; 40000];
        let result = GlobalNtAligner::new(
            NtAlignmentConfig::new(1, -1, -1, -1),
            long_seq
        );
        assert!(matches!(result, Err(AlignmentError::SequenceTooLong)));
    }

    #[test]
    fn test_batch_alignment() {
        let mut al = aligner(b"AGCT");
        let subjects = vec![b"AGCT".as_slice(), b"AGAT".as_slice(), b"".as_slice()];
        let results = al.align_batch(&subjects).unwrap();
        assert_eq!(results.len(), 3);
        let scores: Vec<Score> = results.into_iter().map(|r| r.score).collect();
        assert_eq!(scores, Vec::from([4, 2, -4]));
    }

    #[test]
    fn test_batch_alignment_different_penalties() {
        // match=1, mismatch=-1, gap=-2
        let config = NtAlignmentConfig::new(1, -1, -2, -2);
        let mut al = GlobalNtAligner::new(config, b"A".to_vec()).unwrap();
        
        // Subject 1: "A" (len 1) -> Score: 1 (match)
        // Subject 2: "" (len 0) -> Score: -2 (1 gap in reference)
        // Subject 3: "AA" (len 2) -> Score: 1 - 2 = -1 (match + 1 gap in subject)
        let subjects = vec![b"A".as_slice(), b"".as_slice(), b"AA".as_slice()];
        let results = al.align_batch(&subjects).unwrap();
        let scores: Vec<Score> = results.into_iter().map(|r| r.score).collect();
        
        // If the bug exists:
        // Row 0: [0, -2]
        // Row 1: [ -2, 1] (Sub 1 matches A, Sub 2 treated as null matches nothing)
        // Row 2: [ -4, -1] (Sub 3 matches A at row 1, then gap at row 2)
        // handle_empty_subjects (after loop) reads Row 2, Col 1 for Sub 2 -> -1
        // Correct score for Sub 2 should be Row 0, Col 1 -> -2
        assert_eq!(scores, Vec::from([1, -2, -1]));
    }
}
