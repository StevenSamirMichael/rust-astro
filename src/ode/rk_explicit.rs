use super::types::{ODEState, ODESystem};
use std::fmt::Debug;
use std::ops::{Add, Mul};

pub trait RK<const N: usize> {
    const A: [[f64; N]; N];
    const C: [f64; N];
    const B: [f64; N];

    fn step<F>(x0: f64, y0: &F, h: f64, system: &mut impl ODESystem<F>) -> F
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
        karr.into_iter()
            .enumerate()
            .fold(y0.clone(), |acc, (idx, k)| acc + k * Self::B[idx] * h)
    }

    fn integrate<F>(x0: f64, xend: f64, dx: f64, y0: &F, system: &mut impl ODESystem<F>) -> Vec<F>
    where
        F: ODEState,
    {
        let mut x: f64 = x0;
        let mut v = Vec::<F>::new();
        let mut y = y0.clone();
        while x < xend {
            let ynew = Self::step(x, &y, dx, system);
            v.push(ynew.clone());
            x += dx;
            y = ynew;
        }
        v
    }
}

pub struct RK4 {}
///
/// Buchter tableau for RK4
impl RK<4> for RK4 {
    const A: [[f64; 4]; 4] = [
        [0.0, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.5, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
    ];
    const B: [f64; 4] = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0];
    const C: [f64; 4] = [0.0, 0.5, 0.5, 1.0];
}

pub struct Midpoint {}
impl RK<2> for Midpoint {
    const A: [[f64; 2]; 2] = [[0.0, 0.0], [0.5, 0.0]];
    const B: [f64; 2] = [0.0, 1.0];
    const C: [f64; 2] = [0.0, 0.5];
}

#[cfg(test)]
mod tests {

    use super::super::teststate::TestState;

    use super::*;

    type State = TestState<2>;

    struct HarmonicOscillator {
        k: f64,
    }
    impl HarmonicOscillator {
        fn new(k: f64) -> HarmonicOscillator {
            HarmonicOscillator { k: k }
        }
    }

    impl ODESystem<State> for HarmonicOscillator {
        fn ydot(&mut self, _x: f64, y: &State) -> State {
            State::new(&[y.inner[1], -self.k * y.inner[0]])
        }
    }

    #[test]
    fn testit() {
        let mut system = HarmonicOscillator::new(1.0);
        let y0 = State::new(&[1.0, 0.0]);

        use std::f64::consts::PI;

        // integrating this harmonic oscillator between 0 and 2PI should return to the
        // original state
        let out2 = RK4::integrate(0.0, 2.0 * PI, 0.0001 * 2.0 * PI, &y0, &mut system);
        assert!((out2.last().unwrap().inner[0] - 1.0).abs() < 1.0e-6);
        assert!(out2.last().unwrap().inner[1].abs() < 1.0e-10);
        println!("last = {:?}", out2.last().unwrap());
        println!("len = {}", out2.len());
    }
}
