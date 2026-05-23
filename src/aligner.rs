use crate::alignment::Alignment;
use crate::config::{AlignmentConfig};
use crate::matrix::{Matrix, Idx, AlignmentError};

pub trait Aligner<C>
    where C: AlignmentConfig {

    fn align(&mut self, subject: &[u8]) -> Result<Alignment, AlignmentError>;

    fn align_batch(&mut self, subjects: &[&[u8]]) -> Vec<Result<Alignment, AlignmentError>>;

    fn check_sizes(&self, subject_len: usize, reference_len: usize) -> Result<(), AlignmentError>;
}
