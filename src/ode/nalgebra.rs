//! Implmeent ODEState for NAlgebra OMatrix (owned matrix)

use super::types::ODEState;

impl<const R: usize, const C: usize> ODEState for nalgebra::SMatrix<f64, R, C> {
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
        self.norm() / (self.ode_nelem() as f64).sqrt()
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

    #[inline(always)]
    fn ode_nelem(&self) -> usize {
        self.ncols() * self.nrows()
    }
}
