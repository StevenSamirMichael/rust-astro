use num_traits::Zero;
use std::ops::{Add, Div, Mul, Sub};

use std::fmt::Debug;
use thiserror::Error;

pub type ODEResult<T> = Result<T, Box<dyn std::error::Error + Send + Sync>>;

pub trait ODEState:
    Add<Output = Self>
    + Sub<Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
    + Clone
    + Sized
    + Debug
{
    // Element-wise divisior of self by other
    fn ode_elem_div(&self, other: &Self) -> Self;

    // Element-wise maximum of self with other
    fn ode_elem_max(&self, other: &Self) -> Self;

    // Euclideian norm
    fn ode_norm(&self) -> f64;

    // Element-wise absolute value
    fn ode_abs(&self) -> Self;

    // Add scalar to each element
    fn ode_scalar_add(&self, s: f64) -> Self;

    // Sum of squares
    fn ode_sumsq(&self) -> f64;

    // Number of elements
    fn ode_nelem(&self) -> usize;
}

pub trait ODESystem {
    type Output: ODEState + Zero;
    fn ydot(&mut self, x: f64, y: &Self::Output) -> ODEResult<Self::Output>;
}

#[derive(Debug, Error)]
pub enum ODEError {
    #[error("Stopped at x = {x}.  Reached maximum of {steps} steps.")]
    MaxStepsReached { x: f64, steps: usize },
    #[error("Step Size is Too Small")]
    StepSizeTooSmall,
    #[error("Step error not finite")]
    StepErrorToSmall,
    #[error("Dense output not provided in solution")]
    NoDenseOutputInSolution,
    #[error("Interpolation exceeds solution bounds")]
    InterpExceedsSolutionBounds,
    #[error("Interpolation not implemented for this integrator")]
    InterpNotImplemented,
}

#[derive(Debug, Clone)]
pub struct DenseOutput<S>
where
    S: ODEState,
{
    pub x: Vec<f64>,
    pub h: Vec<f64>,
    pub yprime: Vec<Vec<S>>,
    pub y: Vec<S>,
}

#[derive(Debug, Clone)]
pub struct ODESolution<S>
where
    S: ODEState,
{
    pub nevals: usize,
    pub naccept: usize,
    pub nreject: usize,
    pub x: f64,
    pub y: S,
    pub dense: Option<DenseOutput<S>>,
}

#[derive(Debug, Clone)]
pub struct ODEInterp<S>
where
    S: ODEState,
{
    pub x: Vec<f64>,
    pub y: Vec<S>,
}
