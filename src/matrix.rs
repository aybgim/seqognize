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
    pub scores: Vec<Score>,
    pub ops: Vec<Op>,
    pub rows: usize,
    pub cols: usize,
}

impl Matrix {
    pub fn of(rows: usize, cols: usize) -> Self {
        // We use a sheared layout where diagonal k = row + col is stored contiguously.
        // Index mapping: (row, col) -> (row + col) * rows + row
        // Total diagonals: (rows-1) + (cols-1) = rows + cols - 2. 
        // Max index is (rows+cols-2)*rows + (rows-1).
        // To be safe, we allocate (rows + cols) * rows.
        let size = (rows + cols) * rows;
        Matrix {
            scores: vec![0; size],
            ops: vec![Op::START; size],
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
    fn linear_idx(&self, row: usize, col: usize) -> usize {
        (row + col) * self.rows + row
    }

    #[inline]
    pub fn set(&mut self, idx: Idx, element: Element) {
        let i = self.linear_idx(idx.0, idx.1);
        self.scores[i] = element.score;
        self.ops[i] = element.op;
    }

    #[inline]
    pub fn get(&self, idx: Idx) -> Element {
        let i = self.linear_idx(idx.0, idx.1);
        Element {
            score: self.scores[i],
            op: self.ops[i],
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
            mtx.set((r, c), elements[r][c]);
        }
    }
    mtx
}

pub fn move_back(element: &Element, position: Idx) -> Idx {
    let (row, column) = position;
    match element.op {
        Op::MATCH => (row - 1, column - 1),
        Op::INSERT => (row - 1, column),
        Op::DELETE => (row, column - 1),
        _ => unreachable!()
    }
}
