use seqognize::config::Score;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufReader;

#[derive(Serialize, Deserialize)]
pub struct TestCase {
    pub length: usize,
    pub mutation_rate: u32,
    pub sequence: String,
    pub score: Score,
    pub aligned_sequences: (String, String, String),
}

#[derive(Serialize, Deserialize)]
pub struct TestSuite {
    pub reference: String,
    pub test_cases: Vec<TestCase>,
}


pub fn read_tests() -> TestSuite {
    let path = format!("{}/synth.json", env!("CARGO_MANIFEST_DIR"));
    let file = File::open(path).expect("Failed to open synth.json");
    let reader = BufReader::new(file);

    let test_suite: TestSuite = serde_json::from_reader(reader).expect("Failed to parse synth.json");
    test_suite
}

#[cfg(test)]
mod tests {
    use seqognize::aligner::Aligner;
    use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
    use crate::tests::read_tests;

    #[test]
    fn test_synth() {
        let test_suite = read_tests();
        let reference = test_suite.reference.as_bytes();

        let mut aligner = GlobalNtAligner::new(
            NtAlignmentConfig::new(1, -1, -1, -1),
            reference.to_vec()
        );

        let mutant_sequences: Vec<&[u8]> = test_suite.test_cases.iter()
            .map(|t| t.sequence.as_bytes())
            .collect();

        let results = aligner.align_batch(&mutant_sequences);

        for (i, result) in results.into_iter().enumerate() {
            let test = &test_suite.test_cases[i];
            let alignment = result.expect("Alignment failed during test");
            assert_eq!(test.score, alignment.score);
            assert_eq!(test.aligned_sequences, alignment.aligned_sequences());
        }
    }
}