use super::types::*;

#[derive(Debug)]
pub struct RungaKutta<'a, const N: usize> {
    c: &'a [f64; N],
    b: &'a [f64; N],
    a: &'a [[f64; N]; N],
    dx: f64,
}

impl<'a, const N: usize> RungaKutta<'a, N> {
    fn new(c: &'a [f64; N], b: &'a [f64; N], a: &'a [[f64; N]; N], dx: f64) -> RungaKutta<'a, N> {
        RungaKutta::<'a, N> {
            c: c,
            b: b,
            a: a,
            dx: dx,
        }
    }
}

pub mod euler {
    use super::RungaKutta;
    pub fn new<'a>(dx: f64) -> RungaKutta<'a, 1> {
        const C: [f64; 1] = [1.0];
        const B: [f64; 1] = [1.0];
        const A: [[f64; 1]; 1] = [[1.0]];
        RungaKutta::new(&C, &B, &A, dx)
    }
}

impl<'a, const N: usize> ODESolver for RungaKutta<'a, N> {
    fn integrate(
        &mut self,
        x_start: f64,
        x_end: f64,
        y0: &[f64],
        s: &mut dyn ODESystem,
    ) -> ODEOutput {

        let mut x = x_start;
        let mut y = Vec::from(y0);

        let mut karr = Vec::<&[f64]>::new();
        while x < x_end {
            karr.push(s.ydot(x, y.as_slice()));
            for ik in 1..self.b.len() {
                let mut ya = y.clone();
                for ia in 0..ik {
                    for iv in 0..y.len() {
                        ya[iv] += self.a[ik - 1][ia] * karr[ik - 1][iv] * self.dx;
                    }
                }
                karr.push(s.ydot(x + self.dx + self.c[ik], ya.as_slice()));
            }

            for ik in 0..karr.len() {
                for iv in 0..y.len() {
                    y[iv] += self.b[ik] * karr[ik][iv] * self.dx;
                }
            }

            x += self.dx;
        }
        ODEOutput::new()
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn testit() {
        const C: [f64; 1] = [1.0];
        const B: [f64; 1] = [1.0];
        const A: [[f64; 1]; 1] = [[1.0]];

        let rk = RungaKutta::new(&C, &B, &A);

        println!("rk = {:?}", rk);
    }
}
