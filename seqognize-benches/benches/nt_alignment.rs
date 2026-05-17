#[macro_use]
extern crate criterion;

use criterion::Criterion;
use seqognize::aligner::Aligner;
use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use seqognize_benches::tests::TestSuite;
use std::fs::File;
use std::io::BufReader;

/// Benchmarks the global nucleotide alignment algorithm using a batch of synthetic data.
///
/// This benchmark aligns 100 mutant sequences (ranging from 10 to 5000 bp) against
/// a fixed 1000 bp reference.
///
/// Note: `sample_size` is set to 10 because the total computational work per iteration
/// (100 alignments) is large and would otherwise exceed default timing limits.
fn nt_alignment_benchmark(c: &mut Criterion) {

    let aligner: GlobalNtAligner = GlobalNtAligner {
        config: NtAlignmentConfig {
            match_score: 1,
            mismatch_penalty: -1,
            subject_gap_penalty: -1,
            reference_gap_penalty: -1,
        }
    };
    
    let path = format!("{}/synth.json", env!("CARGO_MANIFEST_DIR"));
    let file = File::open(path).expect("Failed to open synth.json");
    let reader = BufReader::new(file);

    let test_suite: TestSuite = serde_json::from_reader(reader).expect("Failed to parse synth.json");
    let reference = test_suite.reference.as_bytes();
    let mutants: Vec<&[u8]> = test_suite.test_cases.iter()
        .map(|test| test.sequence.as_bytes())
        .collect();

    let mut group = c.benchmark_group("Alignment");
    group.bench_function("NT alignment batch (100 sequences)", |b| {
        b.iter(|| {
            for mutant in &mutants {
                aligner.align(mutant, &reference);
            }
        })
    });

}

criterion_group!(nt_alignment, nt_alignment_benchmark);
criterion_main!(nt_alignment);
