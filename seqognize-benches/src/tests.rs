use serde::{Deserialize, Serialize};
use seqognize::element::Score;

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

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use seqognize::aligner::Aligner;
    use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
    use crate::tests::TestSuite;

    #[test]
    fn test_synth() {
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

        test_suite.test_cases.iter().for_each(|test| {
            let sequence = test.sequence.as_bytes();
            let alignment = aligner.align(&sequence, &reference);
            assert_eq!(test.score, alignment.score);
            assert_eq!(test.aligned_sequences, alignment.aligned_sequences());
        });
    }
}