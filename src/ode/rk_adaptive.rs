use super::rk_adaptive_settings::RKAdaptiveSettings;
use super::types::*;
use num_traits::identities::Zero;

pub trait RKAdaptive<const N: usize> {
    // Butcher Tableau Coefficients
    const A: [[f64; N]; N];
    const C: [f64; N];
    const B: [f64; N];
    const BERR: [f64; N];

    // order
    const ORDER: usize;

    /// First Same as Last
    /// (first compute of next iteration is same as last compute of last iteration)
    const FSAL: bool;

    fn interpolate<S: ODEState>(
        _sol: &ODESolution<S>,
        _xstart: f64,
        _xend: f64,
        _dx: f64,
    ) -> ODEResult<ODEInterp<S>> {
        Err(Box::new(ODEError::InterpNotImplemented))
    }

    /// Convenience function to perform ODE integration
    /// and interpolate from start to finish of integratin
    /// at fixed intervals
    fn integrate_dense<S: ODESystem>(
        x0: f64,
        x_end: f64,
        dx: f64,
        y0: &S::Output,
        system: &mut S,
        settings: &RKAdaptiveSettings,
    ) -> ODEResult<(ODESolution<S::Output>, ODEInterp<S::Output>)> {
        // Make sure dense output is enabled
        let res = match settings.dense_output {
            true => Self::integrate(x0, x_end, y0, system, settings)?,
            false => {
                let mut sc = (*settings).clone();
                sc.dense_output = true;
                Self::integrate(x0, x_end, y0, system, &sc)?
            }
        };
        // Interpolate the result
        let interp = Self::interpolate(&res, x0, x_end, dx)?;
        Ok((res, interp))
    }

    ///
    /// Runga-Kutta integration
    /// with Proportional-Integral controller
    fn integrate<S: ODESystem>(
        x0: f64,
        x_end: f64,
        y0: &S::Output,
        system: &mut S,
        settings: &RKAdaptiveSettings,
    ) -> ODEResult<ODESolution<S::Output>> {
        let mut nevals: usize = 0;
        let mut naccept: usize = 0;
        let mut nreject: usize = 0;
        let mut x = x0.clone();
        let mut y = y0.clone();

        let mut qold: f64 = 1.0e-4;

        // Take guess at initial stepsize
        let mut h = {
            // Adapted from OrdinaryDiffEq.jl
            let sci = (y0.ode_abs() * settings.relerror).ode_scalar_add(settings.abserror);
            let d0 = y0.ode_elem_div(&sci).ode_norm();
            let ydot0 = system.ydot(x0.clone(), &y0)?;
            let d1 = ydot0.ode_elem_div(&sci).ode_norm();
            let h0 = 0.01 * d0 / d1;
            let y1 = y0.clone() + ydot0.clone() * h0;
            let ydot1 = system.ydot(x0 + h0, &y1)?;
            let d2 = (ydot1 - ydot0).ode_elem_div(&sci).ode_norm() / h0;
            let dmax = f64::max(d1, d2);
            let h1 = match dmax < 1e-15 {
                false => (0.01 / f64::max(d1, d2)).powf(1.0 / (1.0 + Self::ORDER as f64)),
                true => f64::max(1e-6, h0 * 1e-3),
            };
            nevals += 2;
            f64::min(100.0 * h0, h1)
        };

        let mut accepted_steps: Option<DenseOutput<S::Output>> = match settings.dense_output {
            false => None,
            true => Some(DenseOutput {
                x: Vec::new(),
                h: Vec::new(),
                yprime: Vec::new(),
                y: Vec::new(),
            }),
        };

        // OK ... lets integrate!
        while x < x_end {
            if x + h > x_end {
                h = x_end - x;
            }

            let mut karr = Vec::new();
            karr.push(system.ydot(x, &y)?);

            // Create the "k"s
            for k in 1..N {
                karr.push(system.ydot(
                    x0 + h * Self::C[k],
                    &(karr.iter().enumerate().fold(y.clone(), |acc, (idx, ki)| {
                        acc + ki.clone() * Self::A[k][idx] * h
                    })),
                )?);
            }

            // Sum the "k"s
            let ynp1 = karr.iter().enumerate().fold(y.clone(), |acc, (idx, k)| {
                acc + k.clone() * Self::B[idx] * h
            });

            // Compute the "error" state by differencing the p and p* orders
            let yerr = karr
                .iter()
                .enumerate()
                .fold(S::Output::zero(), |acc, (idx, k)| {
                    acc + k.clone() * Self::BERR[idx] * h
                });

            // Compute normalized error
            let enorm = {
                let mut ymax = y.ode_abs().ode_elem_max(&ynp1.ode_abs()) * settings.relerror;
                ymax = ymax.ode_scalar_add(settings.abserror);
                let ydiv = yerr.ode_elem_div(&ymax);
                //ydiv.norm() / ((y.ncols() * y.nrows()) as f64).sqrt()
                //(ydiv.map(|x| x.powf(2.0)).sum() / ((y.ncols() * y.nrows()) as f64)).sqrt()
                (ydiv.ode_sumsq() / ydiv.ode_nelem() as f64).sqrt()
            };
            nevals += N;

            if !enorm.is_finite() {
                return Err(Box::new(ODEError::StepErrorToSmall));
            }

            // Run proportional-integral controller on error
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
                match settings.dense_output {
                    true => {
                        let astep = accepted_steps.as_mut().unwrap();
                        astep.x.push(x);
                        astep.h.push(h);
                        astep.yprime.push(karr);
                        astep.y.push(y.clone());
                    }
                    false => {}
                }

                // Adjust step size
                qold = f64::max(enorm, 1.0e-4);
                x += h;
                y = ynp1;
                h = h / q;

                naccept += 1;
                // If dense output, limit step size
            } else {
                nreject += 1;
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
            naccept: naccept,
            nreject: nreject,
            x: x,
            y: y,
            dense: accepted_steps,
        })
    }
}
