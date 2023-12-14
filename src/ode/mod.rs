use num_traits::Zero;
use std::fmt::Debug;
use std::ops::{Add, Mul};

pub trait ODESystem<F>
where
    F: Add<F, Output = F> + Mul<f64, Output = F> + Clone + Debug + Zero,
{
    fn ydot(&mut self, x: f64, y: &F) -> F;
}

pub trait RK<const N: usize> {
    const A: [[f64; N]; N];
    const C: [f64; N];
    const B: [f64; N];

    fn step<F>(x0: f64, y0: &F, h: f64, system: &mut impl ODESystem<F>) -> F
    where
        F: Add<F, Output = F> + Mul<f64, Output = F> + Clone + Debug + Zero,
    {
        let mut karr = Vec::new();
        karr.push(system.ydot(x0, y0));

        // Create the "k"s
        for k in 1..N {
            karr.push(karr.iter().enumerate().fold(F::zero(), |acc, (idx, f)| {
                acc + system.ydot(
                    x0 + h * Self::C[idx],
                    &(y0.clone() + f.clone() * Self::A[k][idx]),
                )
            }))
        }

        // Sum the "k"s
        karr.into_iter()
            .enumerate()
            .fold(y0.clone(), |acc, (idx, k)| acc + k * Self::B[idx])
    }

    fn integrate<F>(x0: f64, xend: f64, dx: f64, y0: &F, system: &mut impl ODESystem<F>) -> Vec<F>
    where
        F: Add<F, Output = F> + Mul<f64, Output = F> + Clone + Debug + Zero,
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

struct RK4 {}

impl RK<4> for RK4 {
    const A: [[f64; 4]; 4] = [
        [0.0, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.5, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
    ];
    const B: [f64; 4] = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0];
    const C: [f64; 4] = [0.0, 0.5, 0.5, 0.0];
}

#[derive(Clone, Debug)]
struct ODEState<const N: usize> {
    inner: [f64; N],
}

impl<const N: usize> ODEState<N> {
    fn zeros() -> ODEState<N> {
        ODEState::<N> { inner: [0.0; N] }
    }

    fn new(input: &[f64]) -> ODEState<N> {
        ODEState::<N> {
            inner: input.try_into().unwrap(),
        }
    }
}

impl<const N: usize> Zero for ODEState<N> {
    fn zero() -> ODEState<N> {
        ODEState::<N> { inner: [0.0; N] }
    }

    fn is_zero(&self) -> bool {
        for ix in 0..self.inner.len() {
            if self.inner[ix].is_zero() == false {
                return false;
            }
        }
        true
    }
}

impl<const N: usize> Add<Self> for ODEState<N> {
    type Output = ODEState<N>;
    fn add(self, other: Self) -> Self {
        Self {
            inner: {
                let mut s = [0.0; N];
                for ix in 0..N {
                    s[ix] = self.inner[ix] + other.inner[ix]
                }
                s
            },
        }
    }
}

impl<const N: usize> Add<f64> for ODEState<N> {
    type Output = ODEState<N>;
    fn add(self, other: f64) -> Self {
        Self {
            inner: {
                let mut s = [0.0; N];
                for ix in 0..N {
                    s[ix] = self.inner[ix] + other;
                }
                s
            },
        }
    }
}

impl<const N: usize> Mul<f64> for ODEState<N> {
    type Output = ODEState<N>;
    fn mul(self, other: f64) -> Self {
        Self {
            inner: {
                let mut s = [0.0; N];
                for ix in 0..N {
                    s[ix] = self.inner[ix] * other;
                }
                s
            },
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    type State = ODEState<2>;

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
        let mut system = HarmonicOscillator::new(0.5);
        let y0 = State::new(&[1.0, 0.0]);

        let out = RK4::step(0.0, &y0, 1.0, &mut system);
        println!("out = {:?}", out);
    }
}
