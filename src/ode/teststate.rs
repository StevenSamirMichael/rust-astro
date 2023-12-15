use super::types::ODEState;
use std::ops::{Add, Mul};

//#[derive(Clone, Debug)]

pub type TestState<const N: usize> = nalgebra::SVector<f64, N>;

/*

pub struct TestState<const N: usize> {
    pub inner: [f64; N],
}

impl<const N: usize> TestState<N> {
    pub fn new(input: &[f64]) -> TestState<N> {
        TestState::<N> {
            inner: input.try_into().unwrap(),
        }
    }
}

impl<const N: usize> Add<Self> for TestState<N> {
    type Output = TestState<N>;
    fn add(self, other: Self) -> Self {
        Self {
            inner: {
                let mut s = [0.0; N];
                for ix in 0..N {
                    s[ix] = self.inner[ix] + other.inner[ix]
                }
                s
            },
        }
    }
}

impl<const N: usize> Add<f64> for TestState<N> {
    type Output = TestState<N>;
    fn add(self, other: f64) -> Self {
        Self {
            inner: {
                let mut s = [0.0; N];
                for ix in 0..N {
                    s[ix] = self.inner[ix] + other;
                }
                s
            },
        }
    }
}

impl<const N: usize> Mul<f64> for TestState<N> {
    type Output = TestState<N>;
    fn mul(self, other: f64) -> Self {
        Self {
            inner: {
                let mut s = [0.0; N];
                for ix in 0..N {
                    s[ix] = self.inner[ix] * other;
                }
                s
            },
        }
    }
}

impl<const N: usize> Iterator for TestState<N> {
    type Item = f64;


}

impl<const N: usize> ODEState for TestState<N> {
    fn infty_norm(&self) -> f64 {
        self.inner.iter().fold(0.0, |acc, x| {
            let xabs = x.clone().abs();
            match acc.partial_cmp(&xabs).unwrap() {
                std::cmp::Ordering::Equal => acc,
                std::cmp::Ordering::Less => xabs,
                std::cmp::Ordering::Greater => acc,
            }
        })
    }
}
*/
