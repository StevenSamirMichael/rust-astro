use super::rk_adaptive::RKAdaptive;

// File below is auto-generated via python script that parses
// data available on web at:
// https://www.sfu.ca/~jverner/RKV98.IIa.Robust.000000351.081209.CoeffsOnlyFLOAT6040

use super::rkv98_nointerp_table as bt;

use super::types::{ODEError, ODEInterp, ODEResult, ODESolution, ODEState};

pub struct RKV98NoInterp {}

const N: usize = 16;

impl RKAdaptive<N, 1> for RKV98NoInterp {
    const ORDER: usize = 9;

    const FSAL: bool = false;

    const B: [f64; N] = bt::B;

    const C: [f64; N] = bt::C;

    const A: [[f64; N]; N] = bt::A;

    const BI: [[f64; 1]; N] = bt::BI;

    const BERR: [f64; N] = {
        let mut berr = [0.0; N];
        let mut ix: usize = 0;
        while ix < N {
            berr[ix] = Self::B[ix] - bt::BHAT[ix];
            ix += 1;
        }
        berr
    };

    /// Interpolate densely calculated solution onto
    /// values that are evenly spaced in "x"
    ///
    fn interpolate<S: ODEState>(
        _sol: &ODESolution<S>,
        _xstart: f64,
        _xend: f64,
        _dx: f64,
    ) -> ODEResult<ODEInterp<S>> {
        Err(Box::new(ODEError::InterpNotImplemented))
    }
}

#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::super::HarmonicOscillator;
    use super::super::RKAdaptive;
    use super::super::RKAdaptiveSettings;
    use super::RKV98NoInterp;

    type State = nalgebra::Vector2<f64>;

    #[test]
    fn testit() -> ODEResult<()> {
        let mut system = HarmonicOscillator::new(1.0);
        let y0 = State::new(1.0, 0.0);

        use std::f64::consts::PI;

        let mut settings = RKAdaptiveSettings::default();
        settings.dense_output = false;
        settings.abserror = 1e-12;
        settings.relerror = 1e-12;

        let res = RKV98NoInterp::integrate(0.0, 2.0 * PI, &y0, &mut system, &settings)?;
        println!("res = {:?}", res);
        let res2 = RKV98NoInterp::integrate(0.0, -2.0 * PI, &y0, &mut system, &settings)?;
        println!("res2 = {:?}", res2);

        assert!((res.y[0] - 1.0).abs() < 1.0e-11);
        assert!(res.y[1].abs() < 1.0e-11);

        Ok(())
    }
}
