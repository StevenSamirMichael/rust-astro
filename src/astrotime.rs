#[derive(PartialEq, PartialOrd, Copy, Clone)]
pub struct AstroTime {
    mjd_tai: f64,
}

use super::datadir;

use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;

pub const UTC1970: f64 = 40587.;
pub const UTC1972: f64 = 41317.;
pub const TAI1972: f64 = UTC1972 + 10. / 86400.;
pub const UTCGPS0: f64 = 44244.; // 1980-01-06
pub const TAIGPS0: f64 = UTCGPS0 + 19. / 86400.;

pub const JD2MJD: f64 = -2400000.5;
pub const MJD2JD: f64 = 2400000.5;

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
    #[derive(Debug)]
    static ref DELTAAT_NEW: Vec<[u32; 2]> = {
        let path = datadir::get()
            .unwrap_or(PathBuf::from("."))
            .join("leap-seconds.list");
        if !path.is_file() {
            panic!("Cannot open leap seconds file");
        }

        let file = match File::open(&path) {
            Err(why) => panic!("Couldn't open {}: {}", path.display(), why),
            Ok(file) => file,
        };

        let mut leapvec = Vec::<[u32; 2]>::new();
        for line in io::BufReader::new(file).lines() {
            match &line {
                Ok(s) => match &s {
                    v if v.chars().next().unwrap() == '#' => (),
                    v => {
                        let split = v.trim().split_whitespace().collect::<Vec<&str>>();
                        if split.len() >= 2 {
                            let a: u32 = split[0].parse().unwrap_or(0);
                            let b: u32 = split[1].parse().unwrap_or(0);
                            if a != 0 && b != 0 {
                                leapvec.push([a, b]);
                            }
                        }
                    }
                },
                Err(_) => (),
            }
        }

        leapvec
    };
}

impl AstroTime {
    #[inline]
    pub fn from_mjd(mjd_tai: f64) -> AstroTime {
        AstroTime { mjd_tai: mjd_tai }
    }

    #[inline]
    pub fn new() -> AstroTime {
        AstroTime { mjd_tai: JD2MJD }
    }
}

#[cfg(test)]
mod tests {
    use crate::astrotime::DELTAAT_NEW;

    #[test]
    fn datadir() {
        println!("deltaat = {:?}", DELTAAT_NEW[1]);
        assert_eq!(1, 1);
    }
}
