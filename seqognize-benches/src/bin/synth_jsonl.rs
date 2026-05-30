//! Synthetic data generator for benchmarking.
//!
//! This tool generates a 1000 bp reference sequence and 5000 randomly mutated
//! sequences with varying lengths (10-5000 bp) and mutation rates (1/10 to 1/1000).
//!
//! The output is written in JSONL format to `seqognize-benches/synth.jsonl`.
//! - 1st line: `{"reference": "..."}`
//! - Following lines: `{"length": ..., "mutation_rate": ..., "sequence": "...", "score": ..., "aligned_sequences": [...]}`
//!
//! # Usage
//! ```powershell
//! cargo run -p seqognize-benches --bin synth --release
//! ```

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use seqognize::aligner::Aligner;
use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use seqognize_parallel::ParallelAligner;
use std::fs::File;
use std::io::{BufWriter, Write};
use seqognize_benches::tests::TestCase;
use serde_json::json;

const BASES: &[u8] = b"ACGT";
const SEED: u64 = 42;
const NUM_TESTS: usize = 100;

fn main() -> std::io::Result<()> {
    let mut rng = StdRng::seed_from_u64(SEED);

    // Generate reference sequence (1000 bp)
    let reference: Vec<u8> = (0..1000).map(|_| BASES[rng.gen_range(0..4)]).collect();
    let reference_str = String::from_utf8(reference.clone()).expect("Invalid reference");

    let base_aligner = GlobalNtAligner::<_>::new(
        NtAlignmentConfig::new(1, -1, -1, -1),
        reference.clone()
    ).expect("Failed to create aligner");
    
    // Use ParallelAligner to leverage all cores
    let mut aligner = ParallelAligner::new(base_aligner, 128);

    println!("Generating {} mutants...", NUM_TESTS);

    // Generate mutant sequences and metadata in memory
    let mutants: Vec<(usize, u32, Vec<u8>)> = (0..NUM_TESTS).map(|_| {
        let length = rng.gen_range(10..=5000);
        let mutation_rate_inv = rng.gen_range(10..=1000);

        let mut sequence = Vec::with_capacity(length);
        for j in 0..length {
            let original_base = if j < reference.len() {
                reference[j]
            } else {
                BASES[rng.gen_range(0..4)]
            };

            if rng.gen_ratio(1, mutation_rate_inv) {
                let mut new_base = BASES[rng.gen_range(0..4)];
                while new_base == original_base {
                    new_base = BASES[rng.gen_range(0..4)];
                }
                sequence.push(new_base);
            } else {
                sequence.push(original_base);
            }
        }
        (length, mutation_rate_inv, sequence)
    }).collect();

    let file = File::create("seqognize-benches/synth.jsonl")?;
    let mut writer = BufWriter::new(file);

    // Write header subdocument
    let header = json!({ "reference": reference_str });
    serde_json::to_writer(&mut writer, &header)?;
    writer.write_all(b"\n")?;

    println!("Aligning and streaming to synth.jsonl...");

    // Stream alignments in parallel
    let seq_slices: Vec<&[u8]> = mutants.iter().map(|m| m.2.as_slice()).collect();
    let results = aligner.align_stream(seq_slices);

    for (i, result) in results.enumerate() {
        let (length, mutation_rate, sequence) = &mutants[i];
        let alignment = result.expect("Alignment failed during synthesis");
        
        let test_case = TestCase {
            length: *length,
            mutation_rate: *mutation_rate,
            sequence: String::from_utf8(sequence.clone()).expect("Invalid sequence"),
            score: alignment.score,
            aligned_sequences: alignment.aligned_sequences(),
        };
        
        serde_json::to_writer(&mut writer, &test_case)?;
        writer.write_all(b"\n")?;
    }

    writer.flush()?;

    println!("Generated seqognize-benches/synth.jsonl with 1 reference and {} mutants.", NUM_TESTS);
    Ok(())
}
