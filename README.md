# Seqognize

**Seqognize** is a sequence alignment library and toolset developed in Rust. This project is a work in progress aimed at demonstrating the use of idiomatic, high-performance Rust in the field of bioinformatics.

The core of the project is a vectorized implementation of the Needleman-Wunsch global alignment algorithm, specifically optimized for nucleotide sequences. By utilizing SIMD (Single Instruction, Multiple Data) through the `wide` crate, it can process multiple alignment scores in parallel, significantly increasing throughput for batch operations.

## Key Features

- **Vectorized Alignment:** Parallel scoring matrix computation using SIMD (AVX2 support with generic fallbacks).
- **Efficient Memory Usage:** Reusable buffers across alignment batches to minimize heap allocations.
- **Nucleotide Optimization:** Specialized scoring and alignment logic for DNA/RNA sequences.

## CLI Usage

The project includes a command-line interface for performing quick alignments between a reference sequence and a subject sequence.

### Building the CLI
```bash
cargo build --release -p seqognize-cli
```

### Running an Alignment
You can run the CLI using `cargo run` or by executing the binary directly.

```bash
# Basic alignment
cargo run -p seqognize-cli -- -r "AGCT" -s "AGAT"

# Alignment with vertical output (useful for visual inspection of gaps/mismatches)
cargo run -p seqognize-cli -- -r "AGCT" -s "AGT" --vertical

# Customizing scores and penalties
cargo run -p seqognize-cli -- -r "AGCT" -s "AGAT" --match 2 --mismatch -2 --sg -3
```

### Options
- `-r, --ref <reference>`: The fixed reference sequence to align against.
- `-s, --sub <subject>`: The subject sequence to be aligned.
- `-m, --match <score>`: Score for a match (default: 1).
- `-x, --mismatch <penalty>`: Penalty for a mismatch (default: -1).
- `--sg <penalty>`: Subject gap opening penalty (default: -1).
- `--rg <penalty>`: Reference gap opening penalty (default: -1).
- `--vertical`: Print the alignment in a vertical format.

---
*Note: This project is currently a work in progress.*
