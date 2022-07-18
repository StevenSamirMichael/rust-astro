//! Is this where to-level documentation goes?
//! apparently not
//! # asdfasdf

/// AstroTime
pub mod astrotime;
pub mod earth_orientation_params;
pub mod gravity;
pub mod itrfcoord;
pub mod jplephem;
pub mod sgp4;
pub mod tle;

//mod sgp4unit;
// pub mod tle;
pub mod iau2000;
pub mod utils;

#[cfg(test)]
mod tests {

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }

    #[test]
    fn tests() {
        #[derive(Debug)]
        struct TS {
            a: f64,
            b: f64,
        }

        fn modts(ts2: &mut TS) {
            ts2.b = 3.0;
        }

        fn mod2(c: &mut f64) {
            *c = 5.0;
        }

        let mut ts: TS = TS { a: 1.0, b: 2.0 };

        println!("ts = {:?}", ts);
        modts(&mut ts);
        println!("ts 2 = {:?}", ts);
        mod2(&mut ts.a);
        println!("ts 3 = {:?}", ts);
    }
}
