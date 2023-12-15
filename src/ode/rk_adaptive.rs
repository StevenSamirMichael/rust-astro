use super::rk_adaptive_settings::RKAdaptiveSettings;
use super::types::*;
use std::fmt::Debug;
use std::ops::{Add, Mul};

pub struct StepOK<F>
where
    F: ODEState,
{
    x: f64,
    y: F,
}

struct StepFail {}

enum StepResult<F>
where
    F: ODEState,
{
    Success(StepOK<F>),
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

    #[inline]
    fn step<F>(x0: f64, y0: &F, h: f64, system: &mut impl ODESystem<F>) -> StepResult<F>
    where
        F: ODEState,
    {
        let mut karr = Vec::new();
        karr.push(system.ydot(x0, y0));

        // Create the "k"s
        for k in 1..N {
            karr.push(system.ydot(
                x0 + h * Self::C[k],
                &(karr.iter().enumerate().fold(y0.clone(), |acc, (idx, ki)| {
                    acc + ki.clone() * Self::A[k][idx] * h
                })),
            ));
        }

        // Sum the "k"s
        let y = karr.iter().enumerate().fold(y0.clone(), |acc, (idx, k)| {
            acc + k.clone() * Self::B[idx] * h
        });
        let ystar = karr
            .into_iter()
            .enumerate()
            .fold(y0.clone(), |acc, (idx, k)| acc + k * Self::BSTAR[idx] * h);

        let yerr = y.clone() + ystar * -1.0;

        StepResult::Success(StepOK { x: x0 + h, y: y })
    }

    fn integrate<F>(
        x0: f64,
        x_end: f64,
        dx: f64,
        y0: &F,
        system: &mut impl ODESystem<F>,
        settings: &RKAdaptiveSettings,
    ) -> ODEResult<ODESolution<F>>
    where
        F: ODEState,
    {
        let mut nevals: usize = 0;
        let mut x = x0.clone();
        let mut y = y0.clone();
        // OK ... lets integrate!

        let mut stepsize: f64 = dx;
        if stepsize == 0.0 {
            stepsize = x_end - x0;
        }
        let mut h = stepsize;

        while x < x_end {
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
                .fold(y0.clone(), |acc, (idx, k)| acc + k * Self::BSTAR[idx] * h);

            let yerr = ynp1.clone() + ynp1star * -1.0;

            nevals += N;
        }

        Ok(ODESolution {
            nevals: nevals,
            y: y0.clone(),
            dense: None,
        })
    }
}
