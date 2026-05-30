#[macro_use]
extern crate criterion;

use criterion::Criterion;
use seqognize::aligner::Aligner;
use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use seqognize::simd_backend::{WideBackend};
use seqognize_parallel::ParallelAligner;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn get_fasta_path() -> String {
    format!("{}/synth.fasta", env!("CARGO_MANIFEST_DIR"))
}

fn read_reference() -> String {
    let file = File::open(get_fasta_path()).expect("Failed to open synth.fasta");
    let mut lines = BufReader::new(file).lines();
    lines.next(); // Skip >reference
    lines.next().expect("Missing reference sequence").expect("Failed to read reference")
}

struct FastaSubjectIterator<R: BufRead> {
    lines: std::io::Lines<R>,
}

impl<R: BufRead> Iterator for FastaSubjectIterator<R> {
    type Item = String;

    fn next(&mut self) -> Option<Self::Item> {
        let _ = self.lines.next()?; // Skip >header
        self.lines.next().map(|l| l.expect("Failed to read sequence"))
    }
}

fn iterate_subjects() -> FastaSubjectIterator<BufReader<File>> {
    let file = File::open(get_fasta_path()).expect("Failed to open synth.fasta");
    let mut lines = BufReader::new(file).lines();
    lines.next(); // skip >reference
    lines.next(); // skip reference sequence
    FastaSubjectIterator { lines }
}

/// Benchmarks the global nucleotide alignment algorithm using streaming FASTA data.
fn nt_alignment_benchmark(c: &mut Criterion) {
    let reference = read_reference();
    
    let mut aligner = GlobalNtAligner::<_, WideBackend>::new(
        NtAlignmentConfig::new(1, -1, -1, -1),
        reference.as_bytes().to_vec()
    ).expect("Failed to create aligner");

    let mut group = c.benchmark_group("Alignment-Streaming-FASTA");
    group.sample_size(10);
    group.bench_function("NT alignment stream (5000 sequences)", |b| {
        b.iter(|| {
            let subjects = iterate_subjects();
            for result in aligner.align_stream(subjects) {
                let _ = criterion::black_box(result);
            }
        })
    });
    group.finish();
}

/// Benchmarks the parallel nucleotide alignment algorithm using streaming FASTA data.
fn parallel_nt_alignment_benchmark(c: &mut Criterion) {
    let reference = read_reference();

    let base_aligner = GlobalNtAligner::<_, WideBackend>::new(
        NtAlignmentConfig::new(1, -1, -1, -1),
        reference.as_bytes().to_vec()
    ).expect("Failed to create aligner");
    
    let mut aligner = ParallelAligner::new(base_aligner, 320);

    let mut group = c.benchmark_group("Alignment-Parallel-Streaming-FASTA");
    group.sample_size(10);
    group.bench_function("Parallel NT alignment stream (5000 sequences)", |b| {
        b.iter(|| {
            let subjects = iterate_subjects();
            for result in aligner.align_stream(subjects) {
                let _ = criterion::black_box(result);
            }
        })
    });
    group.finish();
}

criterion_group!(nt_alignment, nt_alignment_benchmark, parallel_nt_alignment_benchmark);
criterion_main!(nt_alignment);
