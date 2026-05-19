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
    pub config: NtAlignmentConfig
}

impl From<NtAlignmentConfig> for GlobalNtAligner {
    fn from(config: NtAlignmentConfig) -> Self {
        GlobalNtAligner { config }
    }
}

impl Aligner<NtAlignmentConfig> for GlobalNtAligner {
    fn check_sizes(&self, subject_len: usize, reference_len: usize) -> Result<(), AlignmentError> {
        if subject_len > self.config.get_max_subject_size() || reference_len > self.config.get_max_reference_size() {
            return Err(AlignmentError::SequenceTooLong);
        }
        Ok(())
    }

    fn fill_top_row(&self, mtx: &mut Matrix) {
        let ncols = mtx.ncols();
        let mut acc = 0;
        mtx.set((0, 0), Element { score: 0, op: Op::START });
        for col in 1..ncols {
            acc += self.config.get_subject_gap_opening_penalty(col - 1);
            mtx.set((0, col), deletion(acc));
        }
    }

    fn fill_left_column(&self, mtx: &mut Matrix) {
        let nrows = mtx.nrows();
        let mut acc = 0;
        mtx.set((0, 0), Element { score: 0, op: Op::START });
        for row in 1..nrows {
            acc += self.config.get_reference_gap_opening_penalty(row - 1);
            mtx.set((row, 0), insertion(acc));
        }
    }

    fn fill(&self, mtx: &mut Matrix, subject: &[u8], reference: &[u8]) {
        let nrows = mtx.nrows();
        let ncols = mtx.ncols();

        for k in 2..(nrows + ncols - 1) {
            let i_min = 1.max(if k >= ncols { k - ncols + 1 } else { 1 });
            let i_max = (nrows - 1).min(k - 1);

            let mut row = i_min;
            
            // SIMD loop
            while row + 7 <= i_max {
                let mut scores_match_base = [0i16; 8];
                let mut scores_insert = [0i16; 8];
                let mut scores_delete = [0i16; 8];
                let mut sub_scores = [0i16; 8];

                for i in 0..8 {
                    let r = row + i;
                    let c = k - r;
                    let l_idx = r * ncols + c;
                    
                    scores_match_base[i] = mtx.scores[l_idx - ncols - 1];
                    scores_insert[i] = mtx.scores[l_idx - ncols] + self.config.get_reference_gap_opening_penalty(r);
                    scores_delete[i] = mtx.scores[l_idx - 1] + self.config.get_subject_gap_opening_penalty(c);
                    
                    sub_scores[i] = if subject[r - 1] == reference[c - 1] {
                        self.config.match_score
                    } else {
                        self.config.mismatch_penalty
                    };
                }

                let v_scores_match = i16x8::from(scores_match_base) + i16x8::from(sub_scores);
                let v_scores_insert = i16x8::from(scores_insert);
                let v_scores_delete = i16x8::from(scores_delete);

                let v_max_score = v_scores_match.max(v_scores_insert.max(v_scores_delete));
                
                let mask_match = v_max_score.cmp_eq(v_scores_match);
                let mask_insert = v_max_score.cmp_eq(v_scores_insert) & !mask_match;

                let v_ops = mask_match.blend(i16x8::from(Op::MATCH as i16), 
                                mask_insert.blend(i16x8::from(Op::INSERT as i16), 
                                                i16x8::from(Op::DELETE as i16)));

                let final_scores: [i16; 8] = v_max_score.into();
                let final_ops: [i16; 8] = v_ops.into();

                for i in 0..8 {
                    let r = row + i;
                    let c = k - r;
                    let l_idx = r * ncols + c;
                    mtx.scores[l_idx] = final_scores[i];
                    mtx.ops[l_idx] = unsafe { std::mem::transmute(final_ops[i] as u8) };
                }

                row += 8;
            }

            // Scalar Peeling
            for r in row..=i_max {
                let col = k - r;
                let s = subject[r - 1];
                let r_base = reference[col - 1];

                let row_offset = r * ncols;
                let prev_row_offset = (r - 1) * ncols;

                let score_match = mtx.scores[prev_row_offset + col - 1] +
                    self.config.get_substitution_score((r, col), s, r_base);
                let score_insert = mtx.scores[prev_row_offset + col] +
                    self.config.get_reference_gap_opening_penalty(r);
                let score_delete = mtx.scores[row_offset + col - 1] +
                    self.config.get_subject_gap_opening_penalty(col);

                let (score, op) = if score_match >= score_insert && score_match >= score_delete {
                    (score_match, Op::MATCH)
                } else if score_insert >= score_delete {
                    (score_insert, Op::INSERT)
                } else {
                    (score_delete, Op::DELETE)
                };

                mtx.scores[row_offset + col] = score;
                mtx.ops[row_offset + col] = op;
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
        builder.build(mtx.get(end_index).score)
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

    const ALIGNER: GlobalNtAligner = GlobalNtAligner {
        config: NtAlignmentConfig {
            match_score: 1,
            mismatch_penalty: -1,
            subject_gap_penalty: -1,
            reference_gap_penalty: -1,
            max_reference_size: 16383,
            max_subject_size: 16383,
        }
    };

    #[test]
    fn test_fill_top_row() {
        let mut mtx = matrix::of(2, 3);
        ALIGNER.fill_top_row(&mut mtx);
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
    fn test_fill_left_column() {
        let mut mtx = matrix::of(3, 2);
        ALIGNER.fill_left_column(&mut mtx);
        assert_eq!(
            mtx.get((0, 0)),
            Element::default()
        );
        for i in 1..3 {
            assert_eq!(
                mtx.get((i, 0)),
                insertion(-(i as Score))
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
        ALIGNER.fill(&mut mtx, "A".as_bytes(), "A".as_bytes());
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
            ALIGNER.trace_back(&mtx, (1, 1), "A".as_bytes(), "A".as_bytes()),
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
            ALIGNER.trace_back(&mtx, (1, 0), &['A' as u8], &[]),
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
            ALIGNER.trace_back(&mtx, (0, 1), &[], &['A' as u8]),
            Alignment::from("_", "A", -1)
        );
    }

    #[test]
    fn test_match() {
        assert_eq!(
            ALIGNER.align(b"AGCT", b"AGCT").unwrap(),
            Alignment::from("AGCT", "AGCT", 4)
        )
    }

    #[test]
    fn test_mismatch() {
        assert_eq!(
            ALIGNER.align(b"AGAT", b"AGCT").unwrap(),
            Alignment::from("AGAT", "AGCT", 2)
        )
    }

    #[test]
    fn test_insertion() {
        assert_eq!(
            ALIGNER.align(b"AGCT", b"AGT").unwrap(),
            Alignment::from("AGCT", "AG_T", 2)
        )
    }

    #[test]
    fn test_deletion() {
        assert_eq!(
            ALIGNER.align(b"AGT", b"AGCT").unwrap(),
            Alignment::from("AG_T", "AGCT", 2)
        )
    }

    #[test]
    fn test_double_insertion() {
        assert_eq!(
            ALIGNER.align(b"AGCT", b"AT").unwrap(),
            Alignment::from("AGCT", "A__T", 0)
        )
    }

    #[test]
    fn test_double_deletion() {
        assert_eq!(
            ALIGNER.align(b"AT", b"AGCT").unwrap(),
            Alignment::from("A__T", "AGCT", 0)
        )
    }

    #[test]
    fn test_leading_insertion() {
        assert_eq!(
            ALIGNER.align(b"AGCT", b"GCT").unwrap(),
            Alignment::from("AGCT", "_GCT", 2)
        )
    }

    #[test]
    fn test_leading_deletion() {
        assert_eq!(
            ALIGNER.align(b"GCT", b"AGCT").unwrap(),
            Alignment::from("_GCT", "AGCT", 2)
        )
    }

    #[test]
    fn test_trailing_insertion() {
        assert_eq!(
            ALIGNER.align(b"AGCT", b"AGC").unwrap(),
            Alignment::from("AGCT", "AGC_", 2)
        )
    }

    #[test]
    fn test_trailing_deletion() {
        assert_eq!(
            ALIGNER.align(b"AGC", b"AGCT").unwrap(),
            Alignment::from("AGC_", "AGCT", 2)
        )
    }

    #[test]
    fn test_two_insertions() {
        assert_eq!(
            ALIGNER.align(b"AGCT", b"GT").unwrap(),
            Alignment::from("AGCT", "_G_T", 0)
        )
    }

    #[test]
    fn test_two_deletions() {
        assert_eq!(
            ALIGNER.align(b"AC", b"AGCT").unwrap(),
            Alignment::from("A_C_", "AGCT", 0)
        )
    }

    #[test]
    fn test_empty_subject() {
        assert_eq!(
            ALIGNER.align(b"", b"AGCT").unwrap(),
            Alignment::from("____", "AGCT", -4)
        )
    }

    #[test]
    fn test_empty_reference() {
        assert_eq!(
            ALIGNER.align(b"AGCT", b"").unwrap(),
            Alignment::from("AGCT", "____", -4)
        )
    }

    #[test]
    fn test_oversize_subject() {
        let long_seq = vec![b'A'; 40000];
        let result = ALIGNER.align(&long_seq, b"A");
        assert_eq!(result, Err(AlignmentError::SequenceTooLong));
    }

    #[test]
    fn test_oversize_reference() {
        let long_seq = vec![b'A'; 40000];
        let result = ALIGNER.align(b"A", &long_seq);
        assert_eq!(result, Err(AlignmentError::SequenceTooLong));
    }
}
