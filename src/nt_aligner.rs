use crate::config::AlignmentConfig;
use crate::aligner::{Aligner};
use crate::alignment::{Alignment, AlignmentBuilder};
use crate::matrix::{Matrix, Idx, AlignmentError};
use crate::{matrix};
use crate::element::{Score, Element, Op};

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
        }
    }
}

impl Aligner<NtAlignmentConfig> for GlobalNtAligner {
    fn reference(&self) -> &[u8] {
        &self.reference
    }

    fn check_sizes(&self, subject_len: usize, reference_len: usize) -> Result<(), AlignmentError> {
        if subject_len > self.config.get_max_subject_size() || reference_len > self.config.get_max_reference_size() {
            return Err(AlignmentError::SequenceTooLong);
        }
        Ok(())
    }

    fn fill_top_row(&self, mtx: &mut Matrix) {
        let ncols = mtx.ncols();
        mtx.scores[0].copy_from_slice(&self.top_row_scores);
        mtx.ops[0..ncols].copy_from_slice(&self.top_row_ops);
    }

    fn fill_left_column(&self, _mtx: &mut Matrix) {
        // In row recycling, we update column 0 of EVERY row during the fill loop.
        // But for row 0 initialization, we already did it in fill_top_row.
    }

    fn fill(&self, mtx: &mut Matrix, subject: &[u8], reference: &[u8]) {
        let nrows = mtx.nrows();
        let ncols = mtx.ncols();

        for row in 1..nrows {
            let curr = row % 2;
            let prev = (row - 1) % 2;
            let s = subject[row - 1];

            // Initialize first column of the current row
            let acc_ref = mtx.scores[prev][0] + self.config.get_reference_gap_opening_penalty(row - 1);
            mtx.scores[curr][0] = acc_ref;
            mtx.ops[row * ncols] = Op::INSERT;

            for col in 1..ncols {
                let r = reference[col - 1];
                
                let score_match = mtx.scores[prev][col - 1] +
                    self.config.get_substitution_score((row, col), s, r);
                let score_insert = mtx.scores[prev][col] +
                    self.config.get_reference_gap_opening_penalty(row - 1);
                let score_delete = mtx.scores[curr][col - 1] +
                    self.config.get_subject_gap_opening_penalty(col - 1);

                let (score, op) = if score_match >= score_insert && score_match >= score_delete {
                    (score_match, Op::MATCH)
                } else if score_insert >= score_delete {
                    (score_insert, Op::INSERT)
                } else {
                    (score_delete, Op::DELETE)
                };

                mtx.scores[curr][col] = score;
                mtx.ops[row * ncols + col] = op;
            }
        }
    }

    fn end_idx(&self, mtx: &Matrix) -> Idx {
        (mtx.nrows() - 1, mtx.ncols() - 1)
    }

    fn trace_back(&self, mtx: &Matrix, end_index: Idx, subject: &[u8], reference: &[u8]) -> Alignment {
        let mut builder = AlignmentBuilder::new(subject, reference);
        let mut cursor = end_index;
        while cursor != (0, 0) {
            let element = mtx.get(cursor);
            builder.take(element.op, cursor);
            cursor = matrix::move_back(&element, cursor);
        }
        builder.take(Op::START, cursor);
        // build() expects the final score of the entire alignment
        builder.build(mtx.scores[end_index.0 % 2][end_index.1])
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
    use crate::nt_aligner::{GlobalNtAligner, NtAlignmentConfig, deletion, insertion, substitution};
    use crate::aligner::Aligner;
    use crate::matrix;
    use crate::matrix::AlignmentError;
    use crate::alignment::Alignment;
    use crate::element::{Score, Element};

    fn aligner(reference: &[u8]) -> GlobalNtAligner {
        GlobalNtAligner::new(
            NtAlignmentConfig::new(1, -1, -1, -1),
            reference.to_vec()
        )
    }

    #[test]
    fn test_fill_top_row() {
        let mut mtx = matrix::of(2, 3);
        aligner(b"AA").fill_top_row(&mut mtx);
        assert_eq!(
            mtx.get((0, 0)),
            Element::default()
        );
        for i in 1..3 {
            assert_eq!(
                mtx.get((0, i)),
                deletion(-(i as Score))
            );
        }
    }

    #[test]
    fn test_fill_with_match() {
        let mut mtx = matrix::from_elements(
            [
                [Element::default(), deletion(-1)],
                [insertion(-1), substitution(0)]
            ]
        );
        aligner(b"A").fill(&mut mtx, b"A", b"A");
        assert_eq!(
            mtx.get((1, 1)),
            substitution(1)
        );
    }

    #[test]
    fn test_trace_back_snp() {
        let mtx = matrix::from_elements(
            [
                [Element::default(), deletion(-1)],
                [insertion(-1), substitution(1)]
            ]
        );
        assert_eq!(
            aligner(b"A").trace_back(&mtx, (1, 1), b"A", b"A"),
            Alignment::from("A", "A", 1)
        );
    }

    #[test]
    fn test_trace_back_insertion() {
        let mtx = matrix::from_elements(
            [
                [Element::default()],
                [insertion(-1)]
            ]
        );
        assert_eq!(
            aligner(b"").trace_back(&mtx, (1, 0), &[b'A'], &[]),
            Alignment::from("A", "_", -1)
        );
    }

    #[test]
    fn test_trace_back_deletion() {
        let mtx = matrix::from_elements(
            [
                [Element::default(), deletion(-1)]
            ]
        );
        assert_eq!(
            aligner(b"A").trace_back(&mtx, (0, 1), &[], &[b'A']),
            Alignment::from("_", "A", -1)
        );
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
        let result = aligner(&long_seq).align(b"A");
        assert_eq!(result, Err(AlignmentError::SequenceTooLong));
    }
}
