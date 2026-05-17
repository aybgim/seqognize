//! Synthetic data generator for benchmarking.
//!
//! This tool generates a 1000 bp reference sequence and 100 randomly mutated
//! sequences with varying lengths (10-5000 bp) and mutation rates (1/10 to 1/1000).
//!
//! The output is written in FASTA format to `seqognize-benches/synth.fasta`.
//!
//! # Usage
//! ```powershell
//! cargo run -p seqognize-benches --bin synth
//! ```

use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use std::fs::File;
use std::io::{Write, BufWriter};

const BASES: &[u8] = b"ACGT";
const SEED: u64 = 42;

fn main() -> std::io::Result<()> {
    let mut rng = StdRng::seed_from_u64(SEED);
    let file = File::create("seqognize-benches/synth.fasta")?;
    let mut writer = BufWriter::new(file);

    // Generate reference sequence (1000 bp)
    let reference: Vec<u8> = (0..1000)
        .map(|_| BASES[rng.gen_range(0..4)])
        .collect();

    writeln!(writer, ">reference:1000")?;
    writer.write_all(&reference)?;
    writeln!(writer)?;

    // Generate 100 mutated sequences
    for i in 1..=100 {
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

        writeln!(writer, ">{i}|{length}|1/{mutation_rate_inv}")?;
        writer.write_all(&sequence)?;
        writeln!(writer)?;
    }

    println!("Generated seqognize-benches/synth.fasta with 1 reference and 100 mutants.");
    Ok(())
}
