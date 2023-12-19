use nalgebra::{allocator::Allocator, DefaultAllocator, Dim, OMatrix};

use std::fmt::Debug;
use thiserror::Error;

pub type State<R, C> = OMatrix<f64, R, C>;

pub trait ODESystem<R, C>
where
    R: Dim,
    C: Dim,
    DefaultAllocator: Allocator<f64, R, C>,
{
    fn ydot(&mut self, x: f64, y: &State<R, C>) -> State<R, C>;
}

#[derive(Debug, Error)]
pub enum ODEError {
    #[error("Stopped at x = {x}.  Reached maximum of {steps} steps.")]
    MaxStepsReached { x: f64, steps: usize },
    #[error("Step Size is Too Small")]
    StepSizeTooSmall,
    #[error("Step error not finite")]
    StepErrorToSmall,
}

pub type ODEResult<T> = Result<T, ODEError>;

#[derive(Debug, Clone)]
pub struct DenseOutput<R, C>
where
    R: Dim,
    C: Dim,
    DefaultAllocator: Allocator<f64, R, C>,
{
    pub x: Vec<f64>,
    pub y: Vec<State<R, C>>,
}

#[derive(Debug, Clone)]
pub struct ODESolution<R, C>
where
    R: Dim,
    C: Dim,
    DefaultAllocator: Allocator<f64, R, C>,
{
    pub nevals: usize,
    pub y: State<R, C>,
    pub dense: Option<DenseOutput<R, C>>,
}
