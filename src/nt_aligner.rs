use crate::aligner::{Aligner, AlignmentError};
use crate::alignment::{Alignment, AlignmentBuilder, Idx};
use crate::config::{AlignmentConfig, Score, SimdArray, SimdI16, SIMD_ZERO_ARRAY};
use crate::alignment::Op;
use crate::aligner::AlignmentError::SequenceTooLong;

use std::convert::TryFrom;
use wide::CmpEq;

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
    fn get_substitution_score_v(&self, _pos: (usize, usize), subjects: SimdI16, reference: u8) -> SimdI16 {
        let v_ref = SimdI16::from(reference as i16);
        let v_is_match = subjects.cmp_eq(v_ref);
        v_is_match.blend(SimdI16::from(self.match_score), SimdI16::from(self.mismatch_penalty))
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
    pub scores: [Vec<SimdI16>; 2],
    pub ops: Vec<SimdI16>,
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
            scores: [Vec::new(), Vec::new()],
            ops: Vec::new(),
        })
    }
}

impl<C: AlignmentConfig> Aligner<C> for GlobalNtAligner<C> {
    /// Aligns a single subject sequence against the reference.
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

        for chunk_idx in (0..n_subjects).step_by(8) {
            let chunk_end = (chunk_idx + 8).min(n_subjects);
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
    fn prepare_batch_buffers(&mut self, nrows: usize, ncols: usize) {
        if self.scores[0].len() < ncols {
            self.scores[0].resize(ncols, SimdI16::ZERO);
            self.scores[1].resize(ncols, SimdI16::ZERO);
        }
        if self.ops.len() < nrows * ncols {
            self.ops.resize(nrows * ncols, SimdI16::ZERO);
        }
    }

    /// Initializes the first row of the dynamic programming matrix.
    fn fill_first_row(&mut self, ncols: usize) {
        for col in 0..ncols {
            self.scores[0][col] = SimdI16::from(self.top_row_scores[col]);
            self.ops[col] = SimdI16::from(self.top_row_ops[col] as i16);
        }
    }

    /// Executes the vectorized Needleman-Wunsch fill phase for a batch of sequences.
    ///
    /// # Returns
    /// An array containing the final alignment scores for each lane in the SIMD vector.
    fn fill_matrix(&mut self, chunk_subjects: &[&[u8]], nrows: usize, ncols: usize) -> [i16; 8] {
        let mut final_scores = [0i16; 8];

        // Row 0 is initialized but not computed in the SIMD loop. We capture 
        // empty subject scores here before the rolling buffer (size 2) 
        // overwrites Row 0 (e.g., Row 2, Row 4, etc. reuse scores[0]).
        self.handle_empty_subjects(chunk_subjects, &mut final_scores);

        for row in 1..nrows {
            let v_sub_bases = self.gather_subject_bases(chunk_subjects, row);
            let v_ref_gap = SimdI16::from(self.config.get_reference_gap_opening_penalty(row - 1));
            self.compute_first_col(row, v_ref_gap, ncols);
            for col in 1..ncols {
                self.compute_cell_simd(row, col, v_sub_bases, v_ref_gap, ncols);
            }
            self.capture_finished_scores(chunk_subjects, row, &mut final_scores);
        }

        final_scores
    }

    /// Gathers nucleotide bases from the current row for all subjects in the batch into a SIMD vector.
    #[inline(always)]
    fn gather_subject_bases(&self, chunk_subjects: &[&[u8]], row: usize) -> SimdI16 {
        let mut sub_bases = SIMD_ZERO_ARRAY;
        for (i, sub) in chunk_subjects.iter().enumerate() {
            if row <= sub.len() {
                sub_bases[i] = sub[row - 1] as i16;
            }
        }
        SimdI16::from(sub_bases)
    }

    /// Initialize column 0 for this row
    #[inline(always)]
    fn compute_first_col(&mut self, row: usize, v_ref_gap: SimdI16, ncols: usize) {
        let curr = row % 2;
        let prev = (row - 1) % 2;
        self.scores[curr][0] = self.scores[prev][0] + v_ref_gap;
        self.ops[row * ncols] = SimdI16::from(Op::INSERT as i16);
    }

    /// Performs vectorized cell computation for a specific row and column.
    #[inline(always)]
    fn compute_cell_simd(&mut self, row: usize, col: usize, v_sub_bases: SimdI16, v_ref_gap: SimdI16, ncols: usize) {
        let curr = row % 2;
        let prev = (row - 1) % 2;
        let r = self.reference[col - 1];

        // Use the vectorized interface method to support position-specific scores efficiently.
        let v_sub_score = self.config.get_substitution_score_v((row, col), v_sub_bases, r);

        // Recurrence: NW(i, j) = max(diag + substitution, up + gap, left + gap)
        let v_score_match = self.scores[prev][col - 1] + v_sub_score;
        let v_sub_gap = SimdI16::from(self.config.get_subject_gap_opening_penalty(col - 1));

        let v_score_insert = self.scores[prev][col] + v_ref_gap;
        let v_score_delete = self.scores[curr][col - 1] + v_sub_gap;

        let v_max_score = v_score_match.max(v_score_insert).max(v_score_delete);
        
        // Traceback Encoding: Store the winning operation in each lane.
        let mask_match = v_max_score.cmp_eq(v_score_match);
        let mask_insert = v_max_score.cmp_eq(v_score_insert) & !mask_match;

        // Option A: Clean alignment formatting
        let v_ops = mask_match.blend(
            SimdI16::from(Op::MATCH as i16),
            mask_insert.blend(
                SimdI16::from(Op::INSERT as i16),
                SimdI16::from(Op::DELETE as i16),
            ),
        );

        self.scores[curr][col] = v_max_score;
        self.ops[row * ncols + col] = v_ops;
    }

    /// Captures the scores for sequences that have reached their full length at the current row.
    #[inline(always)]
    fn capture_finished_scores(&self, chunk_subjects: &[&[u8]], row: usize, final_scores: &mut [i16; 8]) {
        let curr = row % 2;
        let ref_len = self.reference.len();
        for (i, sub) in chunk_subjects.iter().enumerate() {
            if row == sub.len() {
                let row_scores: SimdArray = self.scores[curr][ref_len].into();
                final_scores[i] = row_scores[i];
            }
        }
    }

    /// Handles sequences with zero length to ensure they get correct gap-only alignment scores.
    #[inline(always)]
    fn handle_empty_subjects(&self, chunk_subjects: &[&[u8]], final_scores: &mut [i16; 8]) {
        let ref_len = self.reference.len();
        for (i, sub) in chunk_subjects.iter().enumerate() {
            if sub.len() == 0 {
                let row0_scores: SimdArray = self.scores[0][ref_len].into();
                final_scores[i] = row0_scores[i];
            }
        }
    }

    /// Performs scalar traceback for each sequence in the batch to reconstruct the alignment paths.
    fn perform_tracebacks(&self, chunk_subjects: &[&[u8]], final_scores: [i16; 8], ncols: usize, all_results: &mut Vec<Alignment>) {
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

    /// Convert SimdI16 to Op
    #[inline(always)]
    fn to_op(&self, i: usize, l_idx: usize) -> Op {
        let ops_simd: SimdI16 = self.ops[l_idx].into();
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
