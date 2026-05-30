//! Synthetic data generator for benchmarking (FASTA format).
//!
//! This tool generates a 1000 bp reference sequence and 5000 randomly mutated
//! sequences in FASTA format for high-throughput benchmarking.
//!
//! Output: `seqognize-benches/synth.fasta`

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::fs::File;
use std::io::{BufWriter, Write};

const BASES: &[u8] = b"ACGT";
const SEED: u64 = 42;
const NUM_TESTS: usize = 5000;

fn main() -> std::io::Result<()> {
    let mut rng = StdRng::seed_from_u64(SEED);

    // Generate reference sequence (1000 bp)
    let reference: Vec<u8> = (0..1000).map(|_| BASES[rng.gen_range(0..4)]).collect();
    let reference_str = String::from_utf8(reference.clone()).expect("Invalid reference");

    let file = File::create("seqognize-benches/synth.fasta")?;
    let mut writer = BufWriter::new(file);

    // Write reference
    writeln!(writer, ">reference")?;
    writeln!(writer, "{}", reference_str)?;

    println!("Generating {} mutants in FASTA format...", NUM_TESTS);

    for i in 0..NUM_TESTS {
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
        
        let seq_str = String::from_utf8(sequence).expect("Invalid sequence");
        writeln!(writer, ">mutant_{}_len_{}_rate_{}", i, length, mutation_rate_inv)?;
        writeln!(writer, "{}", seq_str)?;
    }

    writer.flush()?;

    println!("Generated seqognize-benches/synth.fasta with 1 reference and {} mutants.", NUM_TESTS);
    Ok(())
}
