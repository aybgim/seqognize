use crate::element::Score;

pub trait AlignmentConfig {
    fn get_substitution_score(&self, pos: (usize, usize), s: u8, r: u8) -> Score;
    fn get_subject_gap_opening_penalty(&self, pos: usize) -> Score;
    fn get_reference_gap_opening_penalty(&self, pos: usize) -> Score;

    fn get_max_reference_size(&self) -> usize;

    fn get_max_subject_size(&self) -> usize;
}
