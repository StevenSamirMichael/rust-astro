#[derive(PartialEq, PartialOrd, Copy, Clone)]
pub struct AstroTime {
    mjd_tai: f64,
}

use super::datadir;

use super::earth_orientation_params as eop;

use std::f64::consts::PI;

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

/*
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
*/

pub enum Scale {
    INVALID = -1,
    UTC = 1,
    TT = 2,
    UT1 = 3,
    TAI = 4,
    GPS = 5,
    TDB = 6,
}

lazy_static::lazy_static! {
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

impl std::fmt::Display for AstroTime {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let (year, mon, day, hour, min, sec) = self.to_datetime();

        write!(
            f,
            "{}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
            year,
            mon,
            day,
            hour,
            min,
            sec.floor()
        )
    }
}

impl std::ops::Add<f64> for AstroTime {
    type Output = Self;
    fn add(self, other: f64) -> Self::Output {
        Self {
            mjd_tai: self.mjd_tai + other,
        }
    }
}

impl std::ops::Sub<f64> for AstroTime {
    type Output = Self;
    fn sub(self, other: f64) -> Self::Output {
        Self {
            mjd_tai: self.mjd_tai - other,
        }
    }
}

impl std::ops::Sub for AstroTime {
    type Output = f64;
    fn sub(self, other: AstroTime) -> f64 {
        self.mjd_tai - other.mjd_tai
    }
}

impl AstroTime {
    #[inline]
    pub fn new() -> AstroTime {
        AstroTime { mjd_tai: JD2MJD }
    }

    pub fn from_mjd(val: f64, scale: Scale) -> AstroTime {
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
            Scale::TDB => AstroTime {
                mjd_tai: {
                    let ttc: f64 = (val - (2451545.0 - 2400000.4)) / 36525.0;
                    val - 0.01657 / 86400.0 * (PI / 180.0 * (628.3076 * ttc + 6.2401)).sin()
                        - 32.184 / 86400.0
                },
            },
            Scale::UT1 => AstroTime { mjd_tai: 0.0 },
            Scale::INVALID => AstroTime { mjd_tai: JD2MJD },
        }
    }

    pub fn from_date(year: u32, month: u32, day: u32) -> AstroTime {
        AstroTime::from_mjd(date2mjd_utc(year, month, day) as f64, Scale::UTC)
    }

    pub fn to_date(&self) -> (u32, u32, u32) {
        mjd_utc2date(self.to_mjd(Scale::UTC))
    }

    pub fn to_datetime(&self) -> (u32, u32, u32, u32, u32, f64) {
        let mjd_utc = self.to_mjd(Scale::UTC);
        let (year, month, day) = mjd_utc2date(mjd_utc);
        let fracofday: f64 = mjd_utc - mjd_utc.floor();
        let mut sec: f64 = fracofday * 86400.0;
        let hour: u32 = std::cmp::min((sec / 3600.0).floor() as u32, 23);
        let min: u32 = std::cmp::min((sec as u32 - hour * 3600) / 60 as u32, 59);
        sec = sec - hour as f64 * 3600.0 - min as f64 * 60.0;

        (year, month, day, hour, min, sec)
    }

    pub fn from_datetime(
        year: u32,
        month: u32,
        day: u32,
        hour: u32,
        min: u32,
        sec: f64,
    ) -> AstroTime {
        let mut mjd: f64 = date2mjd_utc(year, month, day) as f64;
        mjd = mjd + (((min + (hour * 60)) * 60) as f64 + sec) / 86400.0;
        AstroTime::from_mjd(mjd, Scale::UTC)
    }

    pub fn to_mjd(&self, ref scale: Scale) -> f64 {
        match scale {
            &Scale::TAI => self.mjd_tai,
            &Scale::GPS => {
                // GPS tracks TAI - 19 seconds
                // after Jan 1 1980, & prior is
                // undefined, but we'll just set it to UTC
                if self.mjd_tai > TAIGPS0 {
                    self.mjd_tai - 19.0 / 86400.0
                } else {
                    self.mjd_tai + mjd_tai2utc_seconds(self.mjd_tai) / 86400.0
                }
            }
            &Scale::TT => self.mjd_tai + 32.184 / 86400.0,
            &Scale::UT1 => {
                // First convert to UTC
                let utc: f64 = self.mjd_tai + mjd_tai2utc_seconds(self.mjd_tai) / 86400.0;

                // Then convert to UT1
                // First earth-orientation parameter is dut1
                // which is (UT1 - UTC) in seconds
                utc + eop::get_from_mjd_utc(utc).unwrap()[0] / 86400.0
            }

            &Scale::UTC => self.mjd_tai + mjd_tai2utc_seconds(self.mjd_tai) / 86400.0,
            &Scale::INVALID => -1.0,
            &Scale::TDB => {
                let tt: f64 = self.mjd_tai + 32.184 / 86400.0;
                let ttc: f64 = (tt - (2451545.0 - 2400000.4)) / 36525.0;
                tt + 0.001657 / 86400.0 * (PI / 180.0 * (628.3076 * ttc + 6.2401)).sin()
            }
        }
    }

    pub fn to_jd(&self, scale: Scale) -> f64 {
        self.to_mjd(scale) + MJD2JD
    }

    pub fn from_jd(jd: f64, scale: Scale) -> AstroTime {
        AstroTime::from_mjd(jd + JD2MJD, scale)
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

fn mjd_utc2date(mjd_utc: f64) -> (u32, u32, u32) {
    // Chapter 15 "Calendars" section 15.11.3 of the
    // Explanatory Suppliment to the Astronomical Almanac
    const Y: i32 = 4716;
    const J: i32 = 1401;
    const M: i32 = 2;
    const N: i32 = 12;
    const R: i32 = 4;
    const P: i32 = 1461;
    const V: i32 = 3;
    const U: i32 = 5;
    const S: i32 = 153;
    const W: i32 = 2;
    const B: i32 = 274277;
    const C: i32 = -38;

    let jd: i32 = (0.5 + mjd_utc + MJD2JD) as i32;
    let mut f: i32 = jd + J;
    f = f + (((4 * jd + B) / 146097) * 3) / 4 + C;
    let e: i32 = R * f + V;
    let g: i32 = (e % P) / R;
    let h: i32 = U * g + W;

    let day: i32 = (h % S) / U + 1;
    let month: i32 = ((h / S + M) % N) + 1;
    let year: i32 = e / P - Y + (N + M - month) / N;

    (year as u32, month as u32, day as u32)
}

fn date2mjd_utc(year: u32, month: u32, day: u32) -> i32 {
    // Chapter 15 "Calendars" section 15.11.3 of the
    // Explanatory Suppliment to the Astronomical Almanac
    // Algorithm 3
    const Y: i32 = 4716;
    const J: i32 = 1401;
    const M: i32 = 2;
    const N: i32 = 12;
    const R: i32 = 4;
    const P: i32 = 1461;
    const Q: i32 = 0;
    const U: i32 = 5;
    const S: i32 = 153;
    const T: i32 = 2;
    const A: i32 = 184;
    const C: i32 = -38;

    let h: i32 = month as i32 - M;
    let g: i32 = year as i32 + Y - (N - h) / N;
    let f: i32 = (h - 1 + N) % N;
    let e: i32 = (P * g + Q) / R + (day as i32) - 1 - J;

    let mut jdn: i32 = e + (S * f + T) / U;
    jdn = jdn - (3 * ((g + A) / 100)) / 4 - C;

    (jdn - 2400001) as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astrotime::DELTAAT_NEW;

    #[test]
    fn datadir() {
        println!("deltaat = {:?}", DELTAAT_NEW[1]);
        assert_eq!(1, 1);
    }

    #[test]
    fn date2mjd_test() {
        // Truth is pulled from leap-seconds file
        // for these examples
        assert_eq!(59219, date2mjd_utc(2021, 1, 5));
        // Try a "leap day"
        assert_eq!(58908, date2mjd_utc(2020, 2, 29));
    }

    #[test]
    fn testdatetime() {
        let tm: AstroTime = AstroTime::from_datetime(2021, 3, 4, 12, 45, 33.0);
        let (year, mon, day, hour, min, sec) = tm.to_datetime();
        assert_eq!(year, 2021);
        assert_eq!(mon, 3);
        assert_eq!(day, 4);
        assert_eq!(hour, 12);
        assert_eq!(min, 45);
        assert!(((sec - 33.0) / 33.0).abs() < 1.0e-5);
    }

    #[test]
    fn add_test() {
        let tm = AstroTime::from_datetime(2021, 3, 4, 11, 20, 41.0);
        let delta: f64 = 0.5;
        let tm2 = tm + delta;
        let (year, mon, day, hour, min, sec) = tm2.to_datetime();
        assert_eq!(year, 2021);
        assert_eq!(mon, 3);
        assert_eq!(day, 4);
        assert_eq!(hour, 23);
        assert_eq!(min, 20);
        assert!(((sec - 41.0) / 41.0).abs() < 1.0e-5);

        let dcalc: f64 = tm2 - tm;
        assert!(((dcalc - delta) / delta).abs() < 1.0e-5);
    }

    #[test]
    fn mjd_utc2date_test() {
        let (year, mon, day) = mjd_utc2date(59219.0);
        assert_eq!(year, 2021);
        assert_eq!(mon, 1);
        assert_eq!(day, 5);

        let (year2, mon2, day2) = mjd_utc2date(58908.0);
        assert_eq!(year2, 2020);
        assert_eq!(mon2, 2);
        assert_eq!(day2, 29);
    }
}
