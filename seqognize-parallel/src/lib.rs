use seqognize::aligner::{Aligner, AlignmentError};
use seqognize::alignment::Alignment;
use seqognize::config::AlignmentConfig;
use seqognize::simd_backend::SimdBackend;
use rayon::prelude::*;
use std::marker::PhantomData;

/// A parallel aligner decorator that uses Rayon to distribute alignment tasks across multiple threads.
pub struct ParallelAligner<C, B, A>
where
    C: AlignmentConfig<B> + Send + Sync,
    B: SimdBackend + Send + Sync,
    A: Aligner<C, B> + Clone + Send + Sync,
{
    base: A,
    chunk_size: usize,
    _phantom: PhantomData<(C, B)>,
}

impl<C, B, A> ParallelAligner<C, B, A>
where
    C: AlignmentConfig<B> + Send + Sync,
    B: SimdBackend + Send + Sync,
    A: Aligner<C, B> + Clone + Send + Sync,
{
    /// Creates a new `ParallelAligner` by wrapping a base aligner.
    ///
    /// # Arguments
    /// * `base` - The base aligner instance to be cloned for each worker thread.
    /// * `chunk_size` - The number of sequences per parallel task.
    pub fn new(base: A, chunk_size: usize) -> Self {
        Self {
            base,
            chunk_size,
            _phantom: PhantomData,
        }
    }
}

impl<C, B, A> Aligner<C, B> for ParallelAligner<C, B, A>
where
    C: AlignmentConfig<B> + Send + Sync,
    B: SimdBackend + Send + Sync,
    A: Aligner<C, B> + Clone + Send + Sync,
{
    fn align(&mut self, subject: &[u8]) -> Result<Alignment, AlignmentError> {
        // For a single alignment, just delegate to the base aligner.
        self.base.align(subject)
    }

    fn align_batch(&mut self, subjects: &[&[u8]]) -> Result<Vec<Alignment>, AlignmentError> {
        subjects
            .par_chunks(self.chunk_size)
            .map_init(
                || self.base.clone(),
                |local_aligner, chunk| local_aligner.align_batch(chunk),
            )
            .collect::<Result<Vec<Vec<_>>, _>>()
            .map(|v| v.into_iter().flatten().collect())
    }
}

pub mod tests {
    #[cfg(test)]
    mod test {
        use super::super::*;
        use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
        use seqognize::simd_backend::WideBackend;

        fn parallel_aligner(reference: &[u8]) -> ParallelAligner<NtAlignmentConfig, WideBackend, GlobalNtAligner<NtAlignmentConfig, WideBackend>> {
            let base = GlobalNtAligner::<_>::new(
                NtAlignmentConfig::new(1, -1, -1, -1),
                reference.to_vec()
            ).unwrap();
            ParallelAligner::new(base, 2) // Small chunk size for testing
        }

        #[test]
        fn test_parallel_batch() {
            let mut al = parallel_aligner(b"AGCT");
            let subjects = vec![b"AGCT".as_slice(), b"AGAT".as_slice(), b"AGT".as_slice(), b"".as_slice()];
            let results = al.align_batch(&subjects).unwrap();
            assert_eq!(results.len(), 4);
            assert_eq!(results[0].score, 4);
            assert_eq!(results[1].score, 2);
            assert_eq!(results[2].score, 2);
            assert_eq!(results[3].score, -4);
        }
    }
}
