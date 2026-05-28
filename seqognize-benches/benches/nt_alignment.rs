#[macro_use]
extern crate criterion;

use criterion::Criterion;
use seqognize::aligner::Aligner;
use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use seqognize::simd_backend::{WideBackend};
use seqognize_parallel::ParallelAligner;
use seqognize_benches::tests::read_tests;

const MULTIPLIER: usize = 50;

/// Benchmarks the global nucleotide alignment algorithm using a large batch of synthetic data.
fn nt_alignment_benchmark(c: &mut Criterion) {
    let test_suite = read_tests();
    let reference = test_suite.reference.as_bytes();
    let mutants_base: Vec<&[u8]> = test_suite.test_cases.iter()
        .map(|test| test.sequence.as_bytes())
        .collect();

    // Multiply the workload to 5000 sequences
    let mut mutants = Vec::with_capacity(mutants_base.len() * MULTIPLIER);
    for _ in 0..MULTIPLIER {
        mutants.extend_from_slice(&mutants_base);
    }

    let mut aligner = GlobalNtAligner::<_, WideBackend>::new(
        NtAlignmentConfig::new(1, -1, -1, -1),
        reference.to_vec()
    ).expect("Failed to create aligner");

    let mut group = c.benchmark_group("Alignment");
    group.sample_size(10);
    group.bench_function("NT alignment batch (5000 sequences)", |b| {
        b.iter(|| {
            let _ = aligner.align_batch(&mutants);
        })
    });
    group.finish();
}

/// Benchmarks the parallel nucleotide alignment algorithm using the same large workload.
fn parallel_nt_alignment_benchmark(c: &mut Criterion) {
    let test_suite = read_tests();
    let reference = test_suite.reference.as_bytes();
    let mutants_base: Vec<&[u8]> = test_suite.test_cases.iter()
        .map(|test| test.sequence.as_bytes())
        .collect();

    // Multiply the workload to 5000 sequences
    let mut mutants = Vec::with_capacity(mutants_base.len() * MULTIPLIER);
    for _ in 0..MULTIPLIER {
        mutants.extend_from_slice(&mutants_base);
    }

    let base_aligner = GlobalNtAligner::<_, WideBackend>::new(
        NtAlignmentConfig::new(1, -1, -1, -1),
        reference.to_vec()
    ).expect("Failed to create aligner");
    
    // Use a multiple of 8 (e.g., 40 * 8 = 320)
    // to ensure efficiency within each parallel task across common SIMD widths.
    let mut aligner = ParallelAligner::new(base_aligner, 320);

    let mut group = c.benchmark_group("Alignment-Parallel");
    group.sample_size(10);
    group.bench_function("Parallel NT alignment batch (5000 sequences)", |b| {
        b.iter(|| {
            let _ = aligner.align_batch(&mutants);
        })
    });
    group.finish();
}

criterion_group!(nt_alignment, nt_alignment_benchmark, parallel_nt_alignment_benchmark);
criterion_main!(nt_alignment);
