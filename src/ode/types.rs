#[derive(Clone, Debug)]
pub struct ODEOutput {
    pub nsteps: usize,
    pub dense_x: Option<Vec<f64>>,
    pub dense_y: Option<Vec<f64>>,
    pub yout: Option<std::vec::Vec<f64>>,
}

impl ODEOutput {
    pub(super) fn new() -> ODEOutput {
        ODEOutput {
            nsteps: 0,
            dense_x: None,
            dense_y: None,
            yout: None,
        }
    }
}

pub trait ODESystem {
    fn ydot(&self, x: f64, y: &[f64]) -> &[f64];
}

pub trait ODESolver {
    fn integrate(
        &mut self,
        x_start: f64,
        x_end: f64,
        y0: &[f64],
        s: &mut dyn ODESystem,
    ) -> ODEOutput;
}
