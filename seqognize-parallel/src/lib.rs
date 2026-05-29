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

    /// Aligns a batch of subject sequences in parallel using Rayon.
    ///
    /// The input batch is divided into contiguous chunks of size `chunk_size`.
    /// Each chunk is processed by a local clone of the base aligner on a separate
    /// worker thread, ensuring high throughput and cache locality.
    ///
    /// # Arguments
    /// * `subjects` - A slice of nucleotide sequences to be aligned.
    ///
    /// # Returns
    /// `Ok(Vec<Alignment>)` containing results for each input sequence, or `Err(AlignmentError)` if any chunk fails.
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

    /// Aligns a stream of subject sequences in parallel.
    ///
    /// This method macro-batches incoming sequences to saturate the Rayon thread pool
    /// while maintaining a streaming interface. It periodically calls `align_batch`
    /// to process these batches in parallel across all available cores.
    ///
    /// # Arguments
    /// * `subjects` - An iterator yielding sequences to be aligned.
    ///
    /// # Returns
    /// An iterator yielding `Result<Alignment, AlignmentError>` for each sequence.
    fn align_stream<'a, I, T>(&'a mut self, subjects: I) -> impl Iterator<Item = Result<Alignment, AlignmentError>> + 'a
    where
        I: IntoIterator<Item = T> + 'a,
        T: AsRef<[u8]> + 'a,
    {
        let mut iterator = subjects.into_iter();
        // Calculate a macro-batch size that is large enough to saturate the thread pool.
        let macro_batch_size = self.chunk_size * rayon::current_num_threads().max(1) * 4;
        std::iter::from_fn(move || {
            let mut batch = Vec::with_capacity(macro_batch_size);
            for _ in 0..macro_batch_size {
                if let Some(s) = iterator.next() {
                    batch.push(s);
                } else {
                    break;
                }
            }

            if batch.is_empty() {
                return None;
            }

            let slices: Vec<&[u8]> = batch.iter().map(|s| s.as_ref()).collect();
            match self.align_batch(&slices) {
                Ok(results) => Some(results.into_iter().map(Ok).collect::<Vec<_>>()),
                Err(e) => Some(vec![Err(e); batch.len()]),
            }
        })
        .flatten()
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

        #[test]
        fn test_parallel_stream() {
            let mut al = parallel_aligner(b"AGCT");
            let subjects = vec!["AGCT", "AGAT", "AGT", ""];
            let results: Vec<_> = al.align_stream(subjects).map(|r| r.unwrap().score).collect();
            assert_eq!(results, vec![4, 2, 2, -4]);
        }
    }
}
