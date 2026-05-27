#[macro_use]
extern crate criterion;

use criterion::Criterion;
use seqognize::aligner::Aligner;
use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use seqognize::simd_backend::WideBackend;
use seqognize_benches::tests::read_tests;

/// Benchmarks the global nucleotide alignment algorithm using a batch of synthetic data.
///
/// This benchmark aligns 100 mutant sequences (ranging from 10 to 5000 bp) against
/// a fixed 1000 bp reference.
///
/// Note: `sample_size` is set to 10 because the total computational work per iteration
/// (100 alignments) is large and would otherwise exceed default timing limits.
fn nt_alignment_benchmark(c: &mut Criterion) {
    let test_suite = read_tests();
    let reference = test_suite.reference.as_bytes();
    let mutants: Vec<&[u8]> = test_suite.test_cases.iter()
        .map(|test| test.sequence.as_bytes())
        .collect();

    let mut aligner = GlobalNtAligner::<_, WideBackend>::new(
        NtAlignmentConfig::new(1, -1, -1, -1),
        reference.to_vec()
    ).expect("Failed to create aligner");

    let mut group = c.benchmark_group("Alignment");
    group.sample_size(10);
    group.bench_function("NT alignment batch (100 sequences)", |b| {
        b.iter(|| {
            let _ = aligner.align_batch(&mutants);
        })
    });
    group.finish();
}

criterion_group!(nt_alignment, nt_alignment_benchmark);
criterion_main!(nt_alignment);
