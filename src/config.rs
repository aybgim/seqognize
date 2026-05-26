
#[cfg(target_feature = "avx2")]
pub use crate::config::avx2_config::*;

#[cfg(target_feature = "avx2")]
mod avx2_config {
    pub type SimdI16 = wide::i16x16;
    pub const LANES: usize = 16;
}

#[cfg(not(target_feature = "avx2"))]
pub use fallback_config::*;

#[cfg(not(target_feature = "avx2"))]
mod fallback_config {
    pub type SimdI16 = wide::i16x8;
    pub const LANES: usize = 8;
}

pub type Score = i16;


pub trait AlignmentConfig {
    fn get_substitution_score(&self, pos: (usize, usize), s: u8, r: u8) -> Score;
    
    #[inline(always)]
    fn get_substitution_score_v(&self, pos: (usize, usize), subjects: SimdI16, reference: u8) -> SimdI16 {
        let subjects_arr: [i16; LANES] = subjects.into();
        let mut results = [0i16; LANES];
        for i in 0..LANES {
            results[i] = self.get_substitution_score(pos, subjects_arr[i] as u8, reference);
        }
        SimdI16::from(results)
    }

    fn get_subject_gap_opening_penalty(&self, pos: usize) -> Score;
    fn get_reference_gap_opening_penalty(&self, pos: usize) -> Score;

    fn get_max_reference_size(&self) -> usize;

    fn get_max_subject_size(&self) -> usize;
}