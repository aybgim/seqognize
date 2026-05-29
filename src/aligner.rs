use std::fmt;
use crate::alignment::Alignment;
use crate::config::AlignmentConfig;
use crate::simd_backend::{SimdBackend, WideBackend};

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum AlignmentError {
    SequenceTooLong,
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentError::SequenceTooLong => write!(f, "Sequences are too long for i16 score range"),
        }
    }
}

pub trait Aligner<C, B = WideBackend>
    where C: AlignmentConfig<B>, B: SimdBackend {

    fn align(&mut self, subject: &[u8]) -> Result<Alignment, AlignmentError>;

    fn align_batch(&mut self, subjects: &[&[u8]]) -> Result<Vec<Alignment>, AlignmentError>;

    fn align_stream<'a, I, T>(&'a mut self, subjects: I) -> impl Iterator<Item = Result<Alignment, AlignmentError>> + 'a
        where I: IntoIterator<Item = T> + 'a, T: AsRef<[u8]> + 'a;
}
