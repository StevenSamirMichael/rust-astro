//! Tsitouras Order 5(4) Runge-Kutta integrator
//!
//! See:
//! https://dblp.org/rec/journals/cma/Tsitouras11
//!
//! Note: in paper, sign is reversed on Bhat[7] ...
//! should be -1.0/66.0, not 1.0/66.0
//!
//! Note also: bhat values are actually error values
//! The nomenclature is confusing
//!

use super::rk_adaptive::RKAdaptive;

use super::types::*;

const A32: f64 = 0.3354806554923570;
const A42: f64 = -6.359448489975075;
const A52: f64 = -11.74888356406283;
const A43: f64 = 4.362295432869581;
const A53: f64 = 7.495539342889836;
const A54: f64 = -0.09249506636175525;
const A62: f64 = -12.92096931784711;
const A63: f64 = 8.159367898576159;
const A64: f64 = -0.07158497328140100;
const A65: f64 = -0.02826905039406838;

const BI11: f64 = -1.0530884977290216;
const BI12: f64 = -1.3299890189751412;
const BI13: f64 = -1.4364028541716351;
const BI14: f64 = 0.7139816917074209;

const BI21: f64 = 0.1017;
const BI22: f64 = -2.1966568338249754;
const BI23: f64 = 1.2949852507374631;

const BI31: f64 = 2.490627285651252793;
const BI32: f64 = -2.38535645472061657;
const BI33: f64 = 1.57803468208092486;

const BI41: f64 = -16.54810288924490272;
const BI42: f64 = -1.21712927295533244;
const BI43: f64 = -0.61620406037800089;

const BI51: f64 = 47.37952196281928122;
const BI52: f64 = -1.203071208372362603;
const BI53: f64 = -0.658047292653547382;

const BI61: f64 = -34.87065786149660974;
const BI62: f64 = -1.2;
const BI63: f64 = -2.0 / 3.0;

const BI71: f64 = 2.5;
const BI72: f64 = -1.0;
const BI73: f64 = -0.6;

pub struct RKTS54 {}

impl RKTS54 {}

impl RKAdaptive<7> for RKTS54 {
    const C: [f64; 7] = [0.0, 0.161, 0.327, 0.9, 0.9800255409045097, 1.0, 1.0];

    const B: [f64; 7] = [
        0.09646076681806523,
        0.01,
        0.4798896504144996,
        1.379008574103742,
        -3.290069515436081,
        2.324710524099774,
        0.0,
    ];

    const BERR: [f64; 7] = [
        0.001780011052226,
        0.000816434459657,
        -0.007880878010262,
        0.144711007173263,
        -0.582357165452555,
        0.458082105929187,
        -1.0 / 66.0,
    ];

    const A: [[f64; 7]; 7] = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [Self::C[1], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [Self::C[2] - A32, A32, 0.0, 0.0, 0.0, 0.0, 0.0],
        [Self::C[3] - A42 - A43, A42, A43, 0.0, 0.0, 0.0, 0.0],
        [Self::C[4] - A52 - A53 - A54, A52, A53, A54, 0.0, 0.0, 0.0],
        [
            Self::C[5] - A62 - A63 - A64 - A65,
            A62,
            A63,
            A64,
            A65,
            0.0,
            0.0,
        ],
        Self::B,
    ];

    const ORDER: usize = 5;

    const FSAL: bool = false;

    fn interpolate<S: ODEState>(
        sol: &ODESolution<S>,
        xstart: f64,
        xend: f64,
        dx: f64,
    ) -> ODEResult<ODEInterp<S>> {
        if sol.dense.is_none() {
            return Err(Box::new(ODEError::NoDenseOutputInSolution));
        }
        if sol.x > xend {
            return Err(Box::new(ODEError::InterpExceedsSolutionBounds));
        }
        let dense = sol.dense.as_ref().unwrap();
        if sol.x < dense.x[0] {
            return Err(Box::new(ODEError::InterpExceedsSolutionBounds));
        }
        let n = ((xend - xstart) / dx).ceil() as usize + 1;
        let mut xarr: Vec<f64> = (0..n).map(|v| v as f64 * dx + xstart).collect();
        if *xarr.last().unwrap() > xend {
            xarr.pop();
            xarr.push(xend);
        }

        //let mut lastidx: usize = 0;
        let yarr: Vec<S> = xarr
            .iter()
            .map(|v| {
                let mut idx = match dense.x.iter().position(|x| *x >= *v) {
                    Some(v) => v,
                    None => dense.x.len(),
                };
                if idx > 0 {
                    idx -= 1;
                }

                // t is fractional distance beween x at idx and idx+1
                let t = (*v - dense.x[idx]) / dense.h[idx];
                let biarr: [f64; 7] = [
                    BI11 * t * (t + BI12) * (t * t + BI13 * t + BI14),
                    BI21 * t * t * (t * t + BI22 * t + BI23),
                    BI31 * t * t * (t * t + BI32 * t + BI33),
                    BI41 * t * t * (t + BI42) * (t + BI43),
                    BI51 * t * t * (t + BI52) * (t + BI53),
                    BI61 * t * t * (t + BI62) * (t + BI63),
                    BI71 * t * t * (t + BI72) * (t + BI73),
                ];
                let yarr = dense.yprime[idx]
                    .iter()
                    .enumerate()
                    .fold(dense.y[idx].clone() / dense.h[idx], |acc, (ix, k)| {
                        acc + k.clone() * biarr[ix]
                    });
                yarr * dense.h[idx]
            })
            .collect();

        Ok(ODEInterp::<S> { x: xarr, y: yarr })
    }
}

#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::super::RKAdaptiveSettings;
    use super::*;
    type State = nalgebra::Vector2<f64>;

    struct HarmonicOscillator {
        k: f64,
    }
    impl HarmonicOscillator {
        fn new(k: f64) -> HarmonicOscillator {
            HarmonicOscillator { k: k }
        }
    }

    impl ODESystem for HarmonicOscillator {
        type Output = nalgebra::Vector2<f64>;
        fn ydot(&mut self, _x: f64, y: &Self::Output) -> ODEResult<Self::Output> {
            Ok(nalgebra::Vector2::<f64>::new(y[1], -self.k * y[0]))
        }
    }

    #[test]
    fn testit() -> ODEResult<()> {
        let mut system = HarmonicOscillator::new(1.0);
        let y0 = State::new(1.0, 0.0);

        use std::f64::consts::PI;

        let mut settings = RKAdaptiveSettings::default();
        settings.dense_output = true;
        settings.abserror = 1e-12;
        settings.relerror = 1e-12;

        let (sol, interp) =
            RKTS54::integrate_dense(0.0, PI / 2.0, PI / 2.0 * 0.05, &y0, &mut system, &settings)?;

        println!("sol evals = {}", sol.nevals);
        interp.x.iter().enumerate().for_each(|(idx, x)| {
            // We know the exact solution for the harmonic oscillator
            let exact = x.cos();
            // Compare with the interpolated result
            let diff = exact - interp.y[idx][0];
            // we set abs and rel error to 1e-10, so lets check!
            assert!(diff.abs() < 1e-12);
        });

        Ok(())
    }
}
