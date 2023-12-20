use super::rk_adaptive_settings::RKAdaptiveSettings;
use super::types::*;

use nalgebra::{allocator::Allocator, DefaultAllocator, Dim};

pub struct StepOK<R, C>
where
    R: Dim,
    C: Dim,
    DefaultAllocator: Allocator<f64, R, C>,
{
    x: f64,
    y: State<R, C>,
}

struct StepFail {}

enum StepResult<R, C>
where
    R: Dim,
    C: Dim,
    DefaultAllocator: Allocator<f64, R, C>,
{
    Success(StepOK<R, C>),
    Failure(StepFail),
}

pub trait RKAdaptive<const N: usize> {
    // Butcher Tableau Coefficients
    const A: [[f64; N]; N];
    const C: [f64; N];
    const B: [f64; N];
    const BSTAR: [f64; N];

    // order
    const ORDER: usize;

    // First Same as Last
    // (first compute of next iteration is same as last compute of last iteration)
    const FSAL: bool;

    fn integrate<R, C>(
        x0: f64,
        x_end: f64,
        y0: &State<R, C>,
        dx: Option<f64>,
        system: &mut impl ODESystem<R, C>,
        settings: &RKAdaptiveSettings,
    ) -> ODEResult<ODESolution<R, C>>
    where
        R: Dim,
        C: Dim,
        DefaultAllocator: Allocator<f64, R, C>,
    {
        let mut nevals: usize = 0;
        let mut x = x0.clone();
        let mut y = y0.clone();

        let mut qold: f64 = 1.0e-4;

        // Take guess at initial stepsize
        let h_init = {
            // Adapted from OrdinaryDiffEq.jl
            let sci = (y0.abs() * settings.relerror).add_scalar(settings.abserror);
            let d0 = y0.component_div(&sci).norm();
            let ydot0 = system.ydot(x0.clone(), &y0);
            let d1 = ydot0.component_div(&sci).norm();
            let h0 = 0.01 * d0 / d1;
            let y1 = y0 + h0 * ydot0.clone();
            let ydot1 = system.ydot(x0 + h0, &y1);
            let d2 = (ydot1 - ydot0).component_div(&sci).norm() / h0;
            let dmax = f64::max(d1, d2);
            let h1 = match dmax < 1e-15 {
                false => (0.01 / f64::max(d1, d2)).powf(1.0 / (1.0 + Self::ORDER as f64)),
                true => f64::max(1e-6, h0 * 1e-3),
            };
            nevals += 2;
            f64::min(100.0 * h0, h1)
        };

        let mut h = match dx {
            None => h_init,
            Some(sdx) => f64::min(h_init, sdx),
        };

        let mut dense_output: Option<DenseOutput<R, C>> = match dx {
            None => None,
            Some(_) => Some(DenseOutput {
                x: Vec::new(),
                y: Vec::new(),
            }),
        };

        // OK ... lets integrate!
        while x < x_end {
            if x + h > x_end {
                h = x_end - x;
            }

            let mut karr = Vec::new();
            karr.push(system.ydot(x, &y));

            // Create the "k"s
            for k in 1..N {
                karr.push(system.ydot(
                    x0 + h * Self::C[k],
                    &(karr.iter().enumerate().fold(y.clone(), |acc, (idx, ki)| {
                        acc + ki.clone() * Self::A[k][idx] * h
                    })),
                ));
            }

            // Sum the "k"s
            let ynp1 = karr.iter().enumerate().fold(y.clone(), |acc, (idx, k)| {
                acc + k.clone() * Self::B[idx] * h
            });
            let ynp1star = karr
                .into_iter()
                .enumerate()
                .fold(y.clone(), |acc, (idx, k)| acc + k * Self::BSTAR[idx] * h);

            let yerr = ynp1.clone() - ynp1star.clone();

            let enorm = {
                let mut ymax = settings.relerror * y0.abs().sup(&ynp1.abs());
                ymax = ymax.add_scalar(settings.abserror);
                let ydiv = yerr.component_div(&ymax);
                //ydiv.norm() / ((y.ncols() * y.nrows()) as f64).sqrt()
                (ydiv.map(|x| x.powf(2.0)).sum() / ((y.ncols() * y.nrows()) as f64)).sqrt()
            };
            nevals += N;

            if !enorm.is_finite() {
                return Err(ODEError::StepErrorToSmall);
            }

            let beta1 = 7.0 / 5.0 / Self::ORDER as f64;
            let beta2 = 2.0 / 5.0 / Self::ORDER as f64;
            let q11 = enorm.powf(beta1);
            let q = {
                let q = q11 / qold.powf(beta2);
                f64::max(
                    1.0 / settings.maxfac,
                    f64::min(1.0 / settings.minfac, q / settings.gamma),
                )
            };

            if (enorm < 1.0) || (h <= settings.dtmin) {
                // If dense output requested, record dense output
                match dx {
                    Some(_) => {
                        let de = dense_output.as_mut().unwrap();
                        de.x.push(x);
                        de.y.push(y);
                    }
                    None => {}
                }

                // Adjust step size
                qold = f64::max(enorm, 1.0e-4);
                x += h;
                y = ynp1;
                h = h / q;

                // If dense output, limit step size
                match dx {
                    Some(sdx) => h = f64::min(h, sdx),
                    None => {}
                }
            } else {
                h = h / f64::min(1.0 / settings.minfac, q11 / settings.gamma);
            }

            /*
            h = h * f64::min(
                settings.maxfac,
                f64::max(
                    settings.minfac,
                    0.9 * (1.0 / enorm).powf(1.0 / (Self::ORDER + 3) as f64),
                ),
            );
            */
        }

        Ok(ODESolution {
            nevals: nevals,
            y: y,
            dense: dense_output,
        })
    }
}
