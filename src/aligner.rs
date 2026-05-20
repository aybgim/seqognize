use crate::alignment::Alignment;
use crate::config::{AlignmentConfig};
use crate::matrix::{Matrix, Idx, AlignmentError};
use crate::matrix;

pub trait Aligner<C>
    where C: AlignmentConfig {

    fn align(&mut self, subject: &[u8]) -> Result<Alignment, AlignmentError>;

    fn align_batch(&mut self, subjects: &[&[u8]]) -> Vec<Result<Alignment, AlignmentError>>;

    fn reference(&self) -> &[u8];

    fn check_sizes(&self, subject_len: usize, reference_len: usize) -> Result<(), AlignmentError>;

    fn fill_top_row(&self, mtx: &mut Matrix);

    fn fill_left_column(&self, mtx: &mut Matrix);

    fn fill(&self, mtx: &mut Matrix, subject: &[u8], reference: &[u8]);

    fn end_idx(&self, mtx: &Matrix) -> Idx;

    fn trace_back(&self, mtx: &Matrix, end_index: Idx, subject: &[u8], reference: &[u8]) -> Alignment;
}
