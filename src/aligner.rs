use std::fmt;
use crate::alignment::Alignment;
use crate::config::AlignmentConfig;

#[derive(Debug, PartialEq)]
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

pub trait Aligner<C>
    where C: AlignmentConfig {

    fn align(&mut self, subject: &[u8]) -> Result<Alignment, AlignmentError>;

    fn align_batch(&mut self, subjects: &[&[u8]]) -> Vec<Result<Alignment, AlignmentError>>;

    fn check_sizes(&self, subject_len: usize, reference_len: usize) -> Result<(), AlignmentError>;
}
