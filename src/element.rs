use std::ops::Add;

pub type Score = i32;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Op {
    START,
    INSERT,
    MATCH,
    DELETE,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Element {
    pub op: Op,
    pub score: Score,
}

impl Add<Score> for Element {
    type Output = Score;

    fn add(self, rhs: Score) -> Self::Output {
        self.score + rhs
    }
}

impl Default for Element {
    fn default() -> Self {
        Element { op: Op::START, score: 0 }
    }
}

