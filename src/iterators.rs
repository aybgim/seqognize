use std::iter::{successors};
use std::ops::Add;
use ndarray::{Dimension};
use ndarray::iter::IterMut;

pub fn accumulate<S, V>(size: usize, supplier: S) -> impl Iterator<Item=V>
    where V: Add<V, Output=V> + Default + Copy,
          S: Fn(usize) -> V {
    let mut range = 0..size;
    successors(
        Some(V::default()),
        move |acc| range.next().map(|n| *acc + supplier(n)),
    )
}

pub fn set_accumulated<V, E>(
    accumulator: impl Iterator<Item=V>,
    setter: IterMut<E, impl Dimension>,
    mapper: fn(V) -> E,
) {
    setter.skip(1)
        .zip(accumulator.skip(1))
        .for_each(|(el, value)| *el = mapper(value));
}