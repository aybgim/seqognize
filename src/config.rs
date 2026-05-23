use crate::element::Score;
use wide::i16x8;

pub trait AlignmentConfig {
    fn get_substitution_score(&self, pos: (usize, usize), s: u8, r: u8) -> Score;
    
    #[inline(always)]
    fn get_substitution_score_v(&self, pos: (usize, usize), subjects: i16x8, reference: u8) -> i16x8 {
        let subjects_arr: [i16; 8] = subjects.into();
        let mut results = [0i16; 8];
        for i in 0..8 {
            results[i] = self.get_substitution_score(pos, subjects_arr[i] as u8, reference);
        }
        i16x8::from(results)
    }

    fn get_subject_gap_opening_penalty(&self, pos: usize) -> Score;
    fn get_reference_gap_opening_penalty(&self, pos: usize) -> Score;

    fn get_max_reference_size(&self) -> usize;

    fn get_max_subject_size(&self) -> usize;
}
