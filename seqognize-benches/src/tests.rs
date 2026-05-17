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