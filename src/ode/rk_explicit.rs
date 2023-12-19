use super::types::{ODESystem, State};
use nalgebra::{allocator::Allocator, DefaultAllocator, Dim};

pub trait RK<const N: usize> {
    const A: [[f64; N]; N];
    const C: [f64; N];
    const B: [f64; N];

    fn step<R, C>(
        x0: f64,
        y0: &State<R, C>,
        h: f64,
        system: &mut impl ODESystem<R, C>,
    ) -> State<R, C>
    where
        R: Dim,
        C: Dim,
        DefaultAllocator: Allocator<f64, R, C>,
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

    fn integrate<R, C>(
        x0: f64,
        xend: f64,
        dx: f64,
        y0: &State<R, C>,
        system: &mut impl ODESystem<R, C>,
    ) -> Vec<State<R, C>>
    where
        R: Dim,
        C: Dim,
        DefaultAllocator: Allocator<f64, R, C>,
    {
        let mut x: f64 = x0;
        let mut v = Vec::new();
        let mut y = y0.clone();
        while x < xend {
            let ynew = Self::step(x, &y, dx, system);
            v.push(ynew.clone());
            x += dx;
            if x > xend {
                x = xend;
            }
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

    impl ODESystem<nalgebra::U2, nalgebra::U1> for HarmonicOscillator {
        fn ydot(&mut self, _x: f64, y: &State) -> State {
            State::new(y[1], -self.k * y[0])
        }
    }

    #[test]
    fn testit() {
        let mut system = HarmonicOscillator::new(1.0);
        let y0 = State::new(1.0, 0.0);

        use std::f64::consts::PI;

        // integrating this harmonic oscillator between 0 and 2PI should return to the
        // original state
        let out2 = RK4::integrate(0.0, 2.0 * PI, 0.0001 * 2.0 * PI, &y0, &mut system);
        assert!((out2.last().unwrap()[0] - 1.0).abs() < 1.0e-6);
        assert!(out2.last().unwrap().abs()[1] < 1.0e-10);
    }
}
