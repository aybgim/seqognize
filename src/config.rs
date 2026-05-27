pub use crate::alignment::Score;
use crate::simd_backend::{SimdBackend, WideBackend};

pub trait AlignmentConfig<B: SimdBackend = WideBackend> {
    fn get_substitution_score(&self, pos: (usize, usize), s: u8, r: u8) -> Score;
    
    #[inline(always)]
    fn get_substitution_score_v(&self, pos: (usize, usize), subjects: B::SimdScore, reference: u8) -> B::SimdScore {
        let subjects_arr = B::vector_to_array(subjects);
        let mut results = B::LanesArray::default();
        for i in 0..B::LANES {
            results[i] = self.get_substitution_score(pos, subjects_arr[i] as u8, reference);
        }
        B::from_array(results)
    }

    fn get_subject_gap_opening_penalty(&self, pos: usize) -> Score;
    fn get_reference_gap_opening_penalty(&self, pos: usize) -> Score;

    fn get_max_reference_size(&self) -> usize;

    fn get_max_subject_size(&self) -> usize;
}
