use std::fmt::Debug;
use std::ops::{Add, Mul};

use thiserror::Error;

pub trait ODEState:
    Add<Output = Self> + Mul<f64, Output = Self> + Clone + Debug + Iterator<Item = f64>
{
}

impl<T> ODEState for T where
    T: Add<Output = T> + Mul<f64, Output = T> + Clone + Debug + Iterator<Item = f64>
{
}

pub trait ODESystem<F>
where
    F: ODEState,
{
    fn ydot(&mut self, x: f64, y: &F) -> F;
}

#[derive(Debug, Error)]
pub enum ODEError {
    #[error("Stopped at x = {x}.  Reached maximum of {steps} steps.")]
    MaxStepsReached { x: f64, steps: usize },
}

pub type ODEResult<T> = Result<T, ODEError>;

pub struct DenseOutput<F>
where
    F: ODEState,
{
    pub x: Vec<f64>,
    pub y: Vec<F>,
}

pub struct ODESolution<F>
where
    F: ODEState,
{
    pub nevals: usize,
    pub y: F,
    pub dense: Option<DenseOutput<F>>,
}
