use super::types::*;

#[derive(Debug)]
pub struct RungaKutta<'a, const N: usize> {
    c: &'a [f64; N],
    b: &'a [f64; N],
    a: &'a [[f64; N]; N],
}

impl<'a, const N: usize> RungaKutta<'a, N> {
    fn new(c: &'a [f64; N], b: &'a [f64; N], a: &'a [[f64; N]; N]) -> RungaKutta<'a, N> {
        RungaKutta::<'a, N> { c: c, b: b, a: a }
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
