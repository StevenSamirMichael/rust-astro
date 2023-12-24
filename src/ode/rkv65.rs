use super::rk_adaptive::RKAdaptive;

use super::types::*;

pub struct RKV65 {}

impl RKAdaptive<10> for RKV65 {
    const ORDER: usize = 6;

    const FSAL: bool = false;

    const C: [f64; 10] = [
        0.0,
        9.0 / 50.0,
        1.0 / 6.0,
        1.0 / 4.0,
        53.0 / 100.0,
        3.0 / 5.0,
        4.0 / 5.0,
        1.0,
        1.0,
        1.0 / 2.0,
    ];

    const B: [f64; 10] = [
        11.0 / 144.0,
        0.0,
        0.0,
        256.0 / 693.0,
        0.0,
        125.0 / 504.0,
        125.0 / 528.0,
        5.0 / 72.0,
        0.0,
        0.0,
    ];

    const BERR: [f64; 10] = {
        const BHAT: [f64; 10] = [
            28.0 / 477.0,
            0.0,
            0.0,
            212.0 / 441.0,
            -312500.0 / 366177.0,
            2125.0 / 1764.0,
            0.0,
            -2105.0 / 35532.0,
            2995.0 / 17766.0,
            0.0,
        ];
        let mut berr = [0.0; 10];
        let mut ix: usize = 0;
        while ix < 6 {
            berr[ix] = Self::B[ix] - BHAT[ix];
            ix += 1;
        }
        berr
    };

    const A: [[f64; 10]; 10] = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [9.0 / 50.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [
            29.0 / 324.0,
            25.0 / 324.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            1.0 / 16.0,
            0.0,
            3.0 / 16.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            79129.0 / 250000.0,
            0.0,
            -261237.0 / 250000.0,
            19663.0 / 15625.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            1336883.0 / 4909125.0,
            0.0,
            -25476.0 / 30875.0,
            194159.0 / 185250.0,
            8225.0 / 78546.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            -2459386.0 / 14727375.0,
            0.0,
            19504.0 / 30875.0,
            2377474.0 / 13615875.0,
            -6157250.0 / 5773131.0,
            902.0 / 735.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            2699.0 / 7410.0,
            0.0,
            -252.0 / 1235.0,
            -1393253.0 / 3993990.0,
            236875.0 / 72618.0,
            -135.0 / 49.0,
            15.0 / 22.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            11.0 / 144.0,
            0.0,
            0.0,
            256.0 / 693.0,
            0.0,
            125.0 / 504.0,
            125.0 / 528.0,
            5.0 / 72.0,
            0.0,
            0.0,
        ],
        [
            1122414499.0 / 13386894336.0,
            0.0,
            0.0,
            1281083.0 / 3946608.0,
            41484375.0 / 162687952.0,
            -184172125.0 / 884040192.0,
            4709875.0 / 84194304.0,
            1019015.0 / 126291456.0,
            -67077.0 / 3508096.0,
            0.0,
        ],
    ];
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
        settings.dense_output = false;
        settings.abserror = 1e-2;
        settings.relerror = 1e-2;

        let sol = RKV65::integrate(0.0, PI / 2.0, &y0, &mut system, &settings)?;

        println!("sol = {:?}", sol);

        Ok(())
    }
}
