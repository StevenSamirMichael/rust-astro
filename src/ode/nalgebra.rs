//! Implmeent ODEState for NAlgebra OMatrix (owned matrix)

use super::types::ODEState;

use nalgebra::{base::allocator::Allocator, DefaultAllocator, DimName};

impl<R: DimName, C: DimName> ODEState for nalgebra::OMatrix<f64, R, C>
where
    DefaultAllocator: Allocator<f64, R, C>,
{
    #[inline(always)]
    fn ode_elem_div(&self, other: &Self) -> Self {
        self.component_div(other)
    }

    #[inline(always)]
    fn ode_elem_max(&self, other: &Self) -> Self {
        self.sup(other)
    }

    #[inline(always)]
    fn ode_norm(&self) -> f64 {
        self.norm()
    }

    #[inline(always)]
    fn ode_sumsq(&self) -> f64 {
        self.map(|x| x * x).sum()
    }

    #[inline(always)]
    fn ode_abs(&self) -> Self {
        self.abs()
    }

    #[inline(always)]
    fn ode_scalar_add(&self, s: f64) -> Self {
        self.add_scalar(s)
    }

    fn ode_nelem(&self) -> usize {
        self.ncols() * self.nrows()
    }
}
