//! Shim for single-threaded rayon replacement
//! based on https://github.com/kornelski/dssim/blob/main/dssim-core/src/lieon.rs
//! Licensed under CC0 and/or AGPLv3 (https://github.com/kornelski/dssim/issues/159)

pub mod prelude {
    pub use super::*;
}

use std::ops::Range;

pub fn current_num_threads() -> usize {
    1
}

pub struct Scope<'scope> {
    _marker: std::marker::PhantomData<&'scope ()>,
}

impl<'scope> Scope<'scope> {
    pub fn spawn<F>(&self, _f: F)
    where
        F: FnOnce(&Scope<'scope>) + Send + 'scope,
    {}
}

pub fn scope<'scope, OP, R>(op: OP) -> R
where
    OP: FnOnce(&Scope<'scope>) -> R + Send,
    R: Send,
{
    let scope = Scope {
        _marker: std::marker::PhantomData,
    };
    op(&scope)
}

pub trait ParIterator: Sized {
    fn with_min_len(self, _one: usize) -> Self { self }
    fn with_max_len(self, _max_len: usize) -> Self { self }
}

impl<T: Iterator> ParIterator for T {
}

pub trait ParSliceLie<T> {
    fn par_chunks_exact(&self, n: usize) -> std::slice::ChunksExact<'_, T>;
    fn par_chunks(&self, n: usize) -> std::slice::Chunks<'_, T>;
}

pub trait ParSliceMutLie<T> {
    fn par_chunks_exact_mut(&mut self, n: usize) -> std::slice::ChunksExactMut<'_, T>;
    fn par_chunks_mut(&mut self, n: usize) -> std::slice::ChunksMut<'_, T>;
}

pub trait ParIntoIterLie<T> {
    type IntoIter;
    fn into_par_iter(self) -> Self::IntoIter;
}

pub trait ParIterLie<'a, T: 'a> {
    type Iter: Iterator<Item = &'a T>;
    fn par_iter(&'a self) -> Self::Iter;
}

pub fn join<A, B>(a: impl FnOnce() -> A, b: impl FnOnce() -> B) -> (A, B) {
    let a = a();
    let b = b();
    (a, b)
}

impl<'a, T> ParSliceLie<T> for &'a [T] {
    fn par_chunks_exact(&self, n: usize) -> std::slice::ChunksExact<'_, T> {
        self.chunks_exact(n)
    }

    fn par_chunks(&self, n: usize) -> std::slice::Chunks<'_, T> {
        self.chunks(n)
    }
}

impl<'a, T> ParSliceLie<T> for &'a mut [T] {
    fn par_chunks_exact(&self, n: usize) -> std::slice::ChunksExact<'_, T> {
        self.chunks_exact(n)
    }

    fn par_chunks(&self, n: usize) -> std::slice::Chunks<'_, T> {
        self.chunks(n)
    }
}

impl<'a, T> ParSliceMutLie<T> for &'a mut [T] {
    fn par_chunks_exact_mut(&mut self, n: usize) -> std::slice::ChunksExactMut<'_, T> {
        self.chunks_exact_mut(n)
    }

    fn par_chunks_mut(&mut self, n: usize) -> std::slice::ChunksMut<'_, T> {
        self.chunks_mut(n)
    }
}

impl<T> ParSliceMutLie<T> for Vec<T> {
    fn par_chunks_exact_mut(&mut self, n: usize) -> std::slice::ChunksExactMut<'_, T> {
        self.as_mut_slice().chunks_exact_mut(n)
    }

    fn par_chunks_mut(&mut self, n: usize) -> std::slice::ChunksMut<'_, T> {
        self.as_mut_slice().chunks_mut(n)
    }
}

impl<'a, T: 'a> ParIterLie<'a, T> for Vec<T> {
    type Iter = std::slice::Iter<'a, T>;

    fn par_iter(&'a self) -> Self::Iter {
        self.as_slice().iter()
    }
}

impl<'a, T: 'a> ParIterLie<'a, T> for &'a [T] {
    type Iter = std::slice::Iter<'a, T>;

    fn par_iter(&'a self) -> Self::Iter {
        self.iter()
    }
}

impl<T> ParIntoIterLie<T> for Vec<T> {
    type IntoIter = std::vec::IntoIter<T>;
    fn into_par_iter(self) -> Self::IntoIter {
        self.into_iter()
    }
}

impl ParIntoIterLie<i32> for Range<i32> {
    type IntoIter = std::ops::Range<i32>;
    fn into_par_iter(self) -> Self::IntoIter {
        self
    }
}

impl ParIntoIterLie<usize> for Range<usize> {
    type IntoIter = std::ops::Range<usize>;
    fn into_par_iter(self) -> Self::IntoIter {
        self
    }
}

pub trait ParIntoIterMutLie<'a, T> {
    type IntoIter;
    fn into_par_iter(self) -> Self::IntoIter;
}

impl<'a, T> ParIntoIterMutLie<'a, T> for &'a mut [T] {
    type IntoIter = std::slice::IterMut<'a, T>;

    fn into_par_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

pub struct ThreadPoolBuilder;

impl ThreadPoolBuilder {
    pub fn new() -> Self {
        ThreadPoolBuilder
    }

    pub fn num_threads(self, _num_cpus: usize) -> Self {
        self
    }

    pub fn build_global(self) -> Result<(), String> {
        Ok(())
    }
}