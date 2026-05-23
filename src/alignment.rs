use core::iter;
use crate::config::Score;

pub const GAP: char = '_';

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Op {
    START,
    INSERT,
    MATCH,
    DELETE,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Idx(pub usize, pub usize);

impl Idx {
    pub fn move_back(&self, op: Op) -> Self {
        match op {
            Op::MATCH => Idx(self.0 - 1, self.1 - 1),
            Op::INSERT => Idx(self.0 - 1, self.1),
            Op::DELETE => Idx(self.0, self.1 - 1),
            _ => unreachable!()
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Anchor {
    pub idx: Idx,
    pub op: Op,
    pub s: u8,
    pub r: u8,
}

impl Anchor {
    const START: Self = Self { idx: Idx(0, 0), op: Op::START, r: 0, s: 0 };

    fn from(idx: Idx, op: Op, s: char, r: char) -> Self {
        Anchor { idx, op, s: s as u8, r: r as u8 }
    }
}

#[derive(Debug, PartialEq)]
pub struct Alignment {
    pub score: Score,
    pub anchors: Vec<Anchor>,
}

impl Alignment {
    pub fn from(subject: &str, reference: &str, score: Score) -> Self {
        Alignment {
            score,
            anchors: to_anchors(subject, reference),
        }
    }

    pub fn pairs(&self, match_symbol: char) -> impl Iterator<Item=(char, char, char)> + '_ {
        self.anchors.iter()
            .rev()
            .skip(1)
            .map(move |a| (
                a.s as char,
                if a.s == a.r { match_symbol } else { ' ' },
                a.r as char
            ))
    }

    pub fn print_horizontal(&self) {
        let als = self.aligned_sequences();
        println!("{}", als.0);
        println!("{}", als.1);
        println!("{}", als.2);
    }

    pub fn print_vertical(&self) {
        self.pairs('-')
            .for_each(|p| println!("{} {} {}", p.0, p.1, p.2));
    }

    pub fn aligned_sequences(&self) -> (String, String, String) {
        let pairs: Vec<(char, char, char)> = self.pairs('|').collect();
        (
            pairs.iter().map(|p| p.0).collect(),
            pairs.iter().map(|p| p.1).collect(),
            pairs.iter().map(|p| p.2).collect()
        )
    }
}

pub struct AlignmentBuilder<'a> {
    anchors: Vec<Anchor>,
    subject: &'a [u8],
    reference: &'a [u8],
}

impl<'a> AlignmentBuilder<'a> {
    pub fn new(subject: &'a [u8], reference: &'a [u8]) -> AlignmentBuilder<'a> {
        AlignmentBuilder {
            anchors: Vec::with_capacity(subject.len() + reference.len()),
            subject,
            reference,
        }
    }

    pub fn take(&mut self, op: Op, idx: Idx) {
        let anchor: Anchor = match op {
            Op::MATCH => Anchor { idx, op, s: self.subject[idx.0 - 1], r: self.reference[idx.1 - 1] },
            Op::DELETE => Anchor { idx, op, s: GAP as u8, r: self.reference[idx.1 - 1] },
            Op::INSERT => Anchor { idx, op, s: self.subject[idx.0 - 1], r: GAP as u8 },
            Op::START => Anchor { idx, op, s: 0, r: 0 }
        };
        self.anchors.push(anchor);
    }

    pub fn build(self, score: Score) -> Alignment {
        Alignment {
            score,
            anchors: self.anchors,
        }
    }
}

fn to_anchors(subject: &str, reference: &str) -> Vec<Anchor> {
    let mut anchors: Vec<Anchor> = iter::once(Anchor::START)
        .chain(from_strings(subject, reference))
        .collect();
    anchors.reverse();
    anchors
}

fn from_strings<'a>(subject: &'a str, reference: &'a str) -> impl Iterator<Item=Anchor> + 'a {
    let mut inc = IdxIncrementer::START;
    subject.chars()
        .zip(reference.chars())
        .map(move |(s, r)|
            Anchor::from(
                inc.with(s, r),
                op(s, r),
                s,
                r,
            )
        )
}

pub fn op(s: char, r: char) -> Op {
    match (s, r) {
        (GAP, _) => Op::DELETE,
        (_, GAP) => Op::INSERT,
        _ => Op::MATCH
    }
}

struct IdxIncrementer {
    s_inc: usize,
    r_inc: usize,
}

impl IdxIncrementer {
    const START: Self = Self { r_inc: 0, s_inc: 0 };

    fn with(&mut self, s: char, r: char) -> Idx {
        Idx(
            Self::with_char(&mut self.s_inc, s),
            Self::with_char(&mut self.r_inc, r)
        )
    }

    fn with_char(i: &mut usize, c: char) -> usize {
        if c != GAP {
            *i += 1;
        }
        *i
    }
}