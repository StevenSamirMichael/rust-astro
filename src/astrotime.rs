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

pub enum Scale {
    INVALID = -1,
    UTC = 1,
    TT = 2,
    UT1 = 3,
    TAI = 4,
    GPS = 5,
    TBD = 6,
}

lazy_static! {
    #[derive(Debug)]
    static ref DELTAAT_NEW: Vec<[u64; 2]> = {
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

        let mut leapvec = Vec::<[u64; 2]>::new();
        for line in io::BufReader::new(file).lines() {
            match &line {
                Ok(s) => match &s {
                    v if v.chars().next().unwrap() == '#' => (),
                    v => {
                        let split = v.trim().split_whitespace().collect::<Vec<&str>>();
                        if split.len() >= 2 {
                            let a: u64 = split[0].parse().unwrap_or(0);
                            let b: u64 = split[1].parse().unwrap_or(0);
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
    pub fn new() -> AstroTime {
        AstroTime { mjd_tai: JD2MJD }
    }

    pub fn from_mjd(val: f64, scale: &Scale) -> AstroTime {
        match scale {
            Scale::TAI => AstroTime {
                mjd_tai: val.clone(),
            },
            Scale::UTC => AstroTime {
                mjd_tai: val + mjd_utc2tai_seconds(val) / 86400.0,
            },
            Scale::TT => AstroTime {
                mjd_tai: val - 32.184 / 86400.0,
            },
            Scale::GPS => AstroTime {
                mjd_tai: {
                    if val > UTCGPS0 {
                        val + 19.0 / 86400.0
                    } else {
                        0.0
                    }
                },
            },
            Scale::TBD => AstroTime { mjd_tai: { 
                let ttc = (val - (2451545.0 - 2400000.4))/36525.0
                0.0 
            } },
            Scale::UT1 => AstroTime { mjd_tai: 0.0 },
            Scale::INVALID => AstroTime { mjd_tai: JD2MJD },
        }
    }
}

/// (TAI - UTC) in seconds for an UTC input modified Julian date
fn mjd_utc2tai_seconds(mjd_utc: f64) -> f64 {
    if mjd_utc > UTC1972 {
        let utc1900: u64 = (mjd_utc as u64 - 15020) * 86400;
        let val = DELTAAT_NEW.iter().find(|&&x| x[0] >= utc1900);
        val.unwrap_or(&[0, 0])[1] as f64
    } else {
        0.0
    }
}

fn mjd_tai2utc_seconds(mjd_tai: f64) -> f64 {
    if mjd_tai > UTC1972 {
        let tai1900: u64 = (mjd_tai as u64 - 15020) * 86400;
        let val = DELTAAT_NEW.iter().find(|&&x| (x[0] + x[1]) > tai1900);
        -(val.unwrap_or(&[0, 0])[1] as f64)
    } else {
        0.0
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
