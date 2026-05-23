#![allow(dead_code)]

use crate::element::{Op, Element, Score};
use std::fmt;

pub type Idx = (usize, usize);

#[derive(Debug, PartialEq)]
pub enum AlignmentError {
    SequenceTooLong,
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentError::SequenceTooLong => write!(f, "Sequences are too long for i16 score range"),
        }
    }
}

pub struct Matrix {
    /// Recycled score rows (current and previous)
    pub scores: [Vec<Score>; 2],
    /// Full operations table for traceback
    pub ops: Vec<Op>,
    pub rows: usize,
    pub cols: usize,
}

impl Matrix {
    pub fn of(rows: usize, cols: usize) -> Self {
        Matrix {
            scores: [vec![0; cols], vec![0; cols]],
            ops: vec![Op::START; rows * cols],
            rows,
            cols,
        }
    }

    pub fn nrows(&self) -> usize {
        self.rows
    }

    pub fn ncols(&self) -> usize {
        self.cols
    }

    #[inline]
    pub fn set_score(&mut self, row_parity: usize, col: usize, score: Score) {
        self.scores[row_parity][col] = score;
    }

    #[inline]
    pub fn get_score(&self, row_parity: usize, col: usize) -> Score {
        self.scores[row_parity][col]
    }

    #[inline]
    pub fn set_op(&mut self, row: usize, col: usize, op: Op) {
        self.ops[row * self.cols + col] = op;
    }

    #[inline]
    pub fn get_op(&self, row: usize, col: usize) -> Op {
        self.ops[row * self.cols + col]
    }

    /// Transparent read for traceback. 
    /// Note: Score is only valid for the very last cell visited (bottom-right).
    #[inline]
    pub fn get(&self, idx: Idx) -> Element {
        Element {
            score: self.scores[idx.0 % 2][idx.1],
            op: self.ops[idx.0 * self.cols + idx.1],
        }
    }
}

pub fn of(num_rows: usize, num_columns: usize) -> Matrix {
    Matrix::of(num_rows, num_columns)
}

pub fn from_elements<const R: usize, const C: usize>(elements: [[Element; C]; R]) -> Matrix {
    let mut mtx = Matrix::of(R, C);
    for r in 0..R {
        for c in 0..C {
            mtx.set_score(r % 2, c, elements[r][c].score);
            mtx.set_op(r, c, elements[r][c].op);
        }
    }
    mtx
}

pub fn move_back_op(op: Op, position: Idx) -> Idx {
    let (row, column) = position;
    match op {
        Op::MATCH => (row - 1, column - 1),
        Op::INSERT => (row - 1, column),
        Op::DELETE => (row, column - 1),
        _ => unreachable!()
    }
}
