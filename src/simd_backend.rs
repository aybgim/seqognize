use crate::alignment::Score;
use wide::CmpEq;
use core::ops::{Index, IndexMut};

/// Trait defining the required SIMD operations for sequence alignment.
pub trait SimdBackend: Copy + Clone {
    /// The SIMD vector type used by this backend.
    type SimdScore: Copy + Clone + Default;
    /// The associated array type for this backend's lanes.
    type LanesArray: Index<usize, Output = Score> + IndexMut<usize> + Default + Copy + Clone;
    /// The number of lanes (SIMD width) for this backend.
    const LANES: usize;

    /// Broadcasts a scalar score to all lanes of a SIMD vector.
    fn splat(score: Score) -> Self::SimdScore;

    /// Performs lane-wise addition of two SIMD vectors.
    fn add(a: Self::SimdScore, b: Self::SimdScore) -> Self::SimdScore;

    /// Performs lane-wise maximum of two SIMD vectors.
    fn max(a: Self::SimdScore, b: Self::SimdScore) -> Self::SimdScore;

    /// Compares two SIMD vectors for equality, returning a mask vector.
    fn cmp_eq(a: Self::SimdScore, b: Self::SimdScore) -> Self::SimdScore;

    /// Performs bitwise AND NOT: (a & !b).
    fn and_not(a: Self::SimdScore, b: Self::SimdScore) -> Self::SimdScore;

    /// Selects lanes from `on_true` where the mask is non-zero, and from `on_false` otherwise.
    fn blend(mask: Self::SimdScore, on_true: Self::SimdScore, on_false: Self::SimdScore) -> Self::SimdScore;

    /// Converts a SIMD vector into an array of scalar scores.
    fn vector_to_array(v: Self::SimdScore) -> Self::LanesArray;

    /// Converts an array of scalar scores into a SIMD vector.
    fn from_array(arr: Self::LanesArray) -> Self::SimdScore;
}

#[cfg(target_feature = "avx2")]
type InternalSimdScore = wide::i16x16;
#[cfg(target_feature = "avx2")]
type InternalLanesArray = [Score; 16];
#[cfg(target_feature = "avx2")]
const INTERNAL_LANES: usize = 16;

#[cfg(not(target_feature = "avx2"))]
type InternalSimdScore = wide::i16x8;
#[cfg(not(target_feature = "avx2"))]
type InternalLanesArray = [Score; 8];
#[cfg(not(target_feature = "avx2"))]
const INTERNAL_LANES: usize = 8;

/// A SIMD backend implementation using the `wide` crate.
#[derive(Copy, Clone, Default)]
pub struct WideBackend;

impl SimdBackend for WideBackend {
    type SimdScore = InternalSimdScore;
    type LanesArray = InternalLanesArray;
    const LANES: usize = INTERNAL_LANES;

    #[inline(always)]
    fn splat(score: Score) -> Self::SimdScore {
        Self::SimdScore::from(score)
    }

    #[inline(always)]
    fn add(a: Self::SimdScore, b: Self::SimdScore) -> Self::SimdScore {
        a + b
    }

    #[inline(always)]
    fn max(a: Self::SimdScore, b: Self::SimdScore) -> Self::SimdScore {
        a.max(b)
    }

    #[inline(always)]
    fn cmp_eq(a: Self::SimdScore, b: Self::SimdScore) -> Self::SimdScore {
        a.cmp_eq(b)
    }

    #[inline(always)]
    fn and_not(a: Self::SimdScore, b: Self::SimdScore) -> Self::SimdScore {
        a & !b
    }

    #[inline(always)]
    fn blend(mask: Self::SimdScore, on_true: Self::SimdScore, on_false: Self::SimdScore) -> Self::SimdScore {
        mask.blend(on_true, on_false)
    }

    #[inline(always)]
    fn vector_to_array(v: Self::SimdScore) -> Self::LanesArray {
        v.into()
    }

    #[inline(always)]
    fn from_array(arr: Self::LanesArray) -> Self::SimdScore {
        Self::SimdScore::from(arr)
    }
}
