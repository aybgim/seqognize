use seqognize::config::Score;
use serde::{Deserialize, Serialize};
use std::fs::File;

#[derive(Serialize, Deserialize)]
pub struct TestCase {
    pub length: usize,
    pub mutation_rate: u32,
    pub sequence: String,
    pub score: Score,
    pub aligned_sequences: (String, String, String),
}

pub fn iterate_tests() -> (String, Box<dyn Iterator<Item = TestCase>>) {
    let path = format!("{}/synth.jsonl", env!("CARGO_MANIFEST_DIR"));
    let file = File::open(path).expect("Failed to open synth.jsonl");
    let mut lines = std::io::BufRead::lines(std::io::BufReader::new(file));

    // Parse header
    let header_line = lines.next().expect("File is empty").expect("Failed to read header");
    let header: serde_json::Value = serde_json::from_str(&header_line).expect("Failed to parse header");
    let reference = header["reference"].as_str().expect("Missing reference in header").to_string();

    // Return reference and an iterator for test cases
    let it = lines.map(|l| {
        let line = l.expect("Failed to read line");
        serde_json::from_str::<TestCase>(&line).expect("Failed to parse test case")
    });

    (reference, Box::new(it))
}

#[cfg(test)]
mod tests {
    use seqognize::aligner::Aligner;
    use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
    use seqognize_parallel::ParallelAligner;
    use seqognize::config::AlignmentConfig;
    use seqognize::simd_backend::SimdBackend;
    use crate::tests::{iterate_tests, TestCase};

    fn test_aligner<C, B, A>(mut aligner: A, test_cases: Box<dyn Iterator<Item = TestCase>>)
    where
        B: SimdBackend,
        C: AlignmentConfig<B>,
        A: Aligner<C, B>,
    {
        let test_cases: Vec<_> = test_cases.collect();
        let seqs: Vec<_> = test_cases.iter().map(|t| t.sequence.as_bytes()).collect();

        for (test, result) in test_cases.iter().zip(aligner.align_stream(seqs)) {
            let alignment = result.expect("Alignment failed");
            assert_eq!(test.score, alignment.score);
            assert_eq!(test.aligned_sequences, alignment.aligned_sequences());
        }
    }

    #[test]
    fn test_synth() {
        let (reference, test_cases) = iterate_tests();

        let aligner = GlobalNtAligner::<_>::new(
            NtAlignmentConfig::new(1, -1, -1, -1),
            reference.as_bytes().to_vec()
        ).expect("Failed to create aligner");

        test_aligner(aligner, test_cases);
    }

    #[test]
    fn test_parallel_synth() {
        let (reference, test_cases) = iterate_tests();

        let base_aligner = GlobalNtAligner::<_>::new(
            NtAlignmentConfig::new(1, -1, -1, -1),
            reference.as_bytes().to_vec()
        ).expect("Failed to create aligner");
        
        let aligner = ParallelAligner::new(base_aligner, 128);

        test_aligner(aligner, test_cases);
    }
}
