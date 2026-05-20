//! Synthetic data generator for benchmarking.
//!
//! This tool generates a 1000 bp reference sequence and 100 randomly mutated
//! sequences with varying lengths (10-5000 bp) and mutation rates (1/10 to 1/1000).
//!
//! The output is written in json format to `seqognize-benches/synth.json`.
//!
//! # Usage
//! ```powershell
//! cargo run -p seqognize-benches --bin synth
//! ```

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use seqognize::aligner::Aligner;
use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use std::fs::File;
use std::io::BufWriter;
use seqognize_benches::tests::{TestCase, TestSuite};

const BASES: &[u8] = b"ACGT";
const SEED: u64 = 42;


const NUM_TESTS: usize = 100;

fn main() -> std::io::Result<()> {
    let mut rng = StdRng::seed_from_u64(SEED);

    // Generate reference sequence (1000 bp)
    let reference: Vec<u8> = (0..1000).map(|_| BASES[rng.gen_range(0..4)]).collect();

    let mut aligner: GlobalNtAligner = GlobalNtAligner::new(
        NtAlignmentConfig::new(1, -1, -1, -1),
        reference.clone()
    );

    let mut test_cases: Vec<TestCase> = Vec::with_capacity(NUM_TESTS);

    // Generate 100 mutated sequences
    for _ in 0..NUM_TESTS {
        let length = rng.gen_range(10..=5000);

        // mutation_rate_inv varies from 10 to 1000, so mutation rate is 1/10 to 1/1000
        let mutation_rate_inv = rng.gen_range(10..=1000);

        let mut sequence = Vec::with_capacity(length);
        for j in 0..length {
            // Pick base from reference if within bounds, otherwise random
            let original_base = if j < reference.len() {
                reference[j]
            } else {
                BASES[rng.gen_range(0..4)]
            };

            // Apply mutation based on rate
            if rng.gen_ratio(1, mutation_rate_inv) {
                // Mutate: pick a different base
                let mut new_base = BASES[rng.gen_range(0..4)];
                while new_base == original_base {
                    new_base = BASES[rng.gen_range(0..4)];
                }
                sequence.push(new_base);
            } else {
                sequence.push(original_base);
            }
        }

        let alignment = aligner.align(&sequence).expect("Alignment failed during synthesis");
        let aligned_sequences = alignment.aligned_sequences();
        test_cases.push(TestCase {
            length,
            mutation_rate: mutation_rate_inv,
            sequence: String::from_utf8(sequence).expect("Invalid reference"),
            score: alignment.score,
            aligned_sequences,
        });
    }

    let file = File::create("seqognize-benches/synth.json")?;
    let writer = BufWriter::new(file);
    let test_suite = TestSuite {
        reference: String::from_utf8(reference).expect("Invalid reference"),
        test_cases,
    };
    serde_json::to_writer_pretty(writer, &test_suite)?;

    println!("Generated seqognize-benches/synth.json with 1 reference and 100 mutants.");
    Ok(())
}
