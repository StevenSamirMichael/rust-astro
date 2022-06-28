#[derive(PartialEq, PartialOrd, Copy, Clone)]
pub struct AstroTime {
    mjd_tai: f64,
}

use super::datadir;

use std::path::PathBuf;

pub const UTC1970: f64 = 40587.;
pub const UTC1972: f64 = 41317.;
pub const TAI1972: f64 = UTC1972 + 10. / 86400.;
pub const UTCGPS0: f64 = 44244.; // 1980-01-06
pub const TAIGPS0: f64 = UTCGPS0 + 19. / 86400.;

const DELTAAT_OLD: [[f64; 4]; 15] = [
    [36204., 0., 36204., 0.],
    [36934., 1.4178180, 37300., 0.001296],
    [37300., 1.4228180, 37300., 0.001296],
    [37512., 1.3728180, 37300., 0.001296],
    [37665., 1.8458580, 37665., 0.0011232],
    [38334., 1.9458580, 37665., 0.0011232],
    [38395., 3.2401300, 38761., 0.001296],
    [38486., 3.3401300, 38761., 0.001296],
    [38639., 3.4401300, 38761., 0.001296],
    [38761., 3.5401300, 38761., 0.001296],
    [38820., 3.6401300, 38761., 0.001296],
    [38942., 3.7401300, 38761., 0.001296],
    [39004., 3.8401300, 38761., 0.001296],
    [39126., 4.3131700, 39126., 0.002592],
    [39887., 4.2131700, 39126., 0.002592],
];

lazy_static! {
    static ref DELTAAT_NEW: Vec<[f64; 4]> = {
        let file = datadir::get()
            .unwrap_or(PathBuf::from("."))
            .join("leap-seconds.list");
        if !file.is_file() {
            panic!("Cannot open leap seconds file");
        }

        let v = Vec::<[f64; 4]>::new();

        v
    };
}

impl AstroTime {
    #[inline]
    pub fn new(mjd_tai: f64) -> AstroTime {
        AstroTime { mjd_tai: mjd_tai }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn datadir() {
        assert_eq!(1, 1);
    }
}
