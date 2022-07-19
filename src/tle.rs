use super::astrotime::AstroTime;
use crate::sgp4::SatRec;

#[derive(PartialEq, PartialOrd, Clone, Debug)]
pub struct TLE {
    /// Name of satellite
    pub name: String,
    /// String describing launch
    pub intl_desig: String,
    /// Satellite NORAD number
    pub sat_num: i32,
    /// Launch year
    pub desig_year: i32,
    /// Numbered launch of year
    pub desig_launch: i32,
    /// Piece of launch
    pub desig_piece: String,
    /// TLE epoch
    pub epoch: AstroTime,
    /// One half of 1st derivative of mean motion wrt time, in revs/day^2
    pub mean_motion_dot: f64,
    /// One sixth of 2nd derivative of mean motion wrt tim, in revs/day^3
    pub mean_motion_dot_dot: f64,
    /// Starred ballistic coefficient, in units of inverse Earth radii
    pub bstar: f64,
    /// Usually 0
    pub ephem_type: u8,
    /// Bulliten number
    pub element_num: i32,
    /// Inclination, degrees
    pub inclination: f64,
    /// Right ascension of ascending node, degrees
    pub raan: f64,
    /// Eccentricity
    pub eccen: f64,
    /// Argument of perigee, degrees
    pub arg_of_perigee: f64,
    /// Mean anomaly, degrees
    pub mean_anomaly: f64,
    /// Mean motion, revs / day
    pub mean_motion: f64,
    /// Revolution number
    pub rev_num: i32,

    pub(in crate) satrec: Option<SatRec>,
}

impl TLE {
    pub fn new() -> TLE {
        TLE {
            name: "none".to_string(),
            intl_desig: "".to_string(),
            sat_num: 0,
            desig_year: 0,
            desig_launch: 0,
            desig_piece: "A".to_string(),
            epoch: AstroTime::new(),
            mean_motion_dot: 0.0,
            mean_motion_dot_dot: 0.0,
            bstar: 0.0,
            ephem_type: 'U' as u8,
            element_num: 0,
            inclination: 0.0,
            raan: 0.0,
            eccen: 0.0,
            arg_of_perigee: 0.0,
            mean_anomaly: 0.0,
            mean_motion: 0.0,
            rev_num: 0,
            satrec: None,
        }
    }

    /// Load 3 lines as string into a structure representing
    /// a Two-Line Element Set  (TLE)
    ///
    /// The TLE can then be used to compute satellite position and
    /// velocity as a function of time.
    ///
    /// For details, see [here](https://en.wikipedia.org/wiki/Two-line_element_set)
    ///
    /// # Arguments:
    ///
    ///   * line0: the "0"th line of the TLE, which sometimes contains
    ///     the satellite name
    ///   * line1: the 1st line of TLE
    ///   * line2: the 2nd line of the TLE
    ///
    /// # Returns:
    ///
    ///  * A TLE object or string indicating error condition
    ///
    /// # Example
    ///
    /// ```
    ///
    /// use astro::TLE;
    /// let line0: &str = "0 INTELSAT 902";
    /// let line1: &str = "1 26900U 01039A   06106.74503247  .00000045  00000-0  10000-3 0  8290";
    /// let line2: &str = "2 26900   0.0164 266.5378 0003319  86.1794 182.2590  1.00273847 16981   9300.";
    /// let tle = TLE::load_3line(&line0.to_string(),
    ///     &line1.to_string(),
    ///     &line2.to_string()
    ///     ).unwrap();
    ///
    /// ```
    ///
    pub fn load_3line(line0: &String, line1: &String, line2: &String) -> Result<TLE, String> {
        match TLE::load_2line(line1, line2) {
            Ok(mut tle) => {
                tle.name = {
                    if line0.len() > 2 && line0.chars().nth(0).unwrap() == '0' {
                        line0[2..].to_string()
                    } else {
                        line0.clone()
                    }
                };
                Ok(tle)
            }
            Err(e) => Err(e),
        }
    }

    pub fn load_2line(line1: &String, line2: &String) -> Result<TLE, String> {
        let mut year: u32 = {
            let mut mstr: String = "1".to_owned();
            mstr.push_str(&line1[18..20]);
            let mut s: u32 = match mstr.parse() {
                Ok(y) => y,
                Err(_) => return Err("Could not parse year".to_string()),
            };
            s = s - 100;
            s
        };
        if year > 57 {
            year += 1900;
        } else {
            year += 2000;
        }
        let day_of_year: f64 = match line1[20..32].parse() {
            Ok(y) => y,
            Err(_) => return Err("Could not parse day of year".to_string()),
        };

        // Note: day_of_year starts from 1, not zero,
        // also, go from Jan 2 to avoid leap-second
        // issues, hence the "-2" at end
        let epoch = AstroTime::from_date(year, 1, 2).add_utc_days(day_of_year - 2.0);

        Ok(TLE {
            name: "none".to_string(),
            sat_num: {
                match line1[2..7].parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse sat number".to_string()),
                }
            },
            intl_desig: { line1[9..16].trim().to_string() },
            desig_year: {
                match line1[9..11].trim().parse() {
                    Ok(l) => l,
                    Err(_) => 70,
                }
            },
            desig_launch: {
                match line1[11..14].trim().parse() {
                    Ok(l) => l,
                    Err(_) => 0,
                }
            },
            desig_piece: {
                match line1[14..18].trim().parse() {
                    Ok(l) => l,
                    Err(_) => return Err("Could not parse desig_piece".to_string()),
                }
            },
            epoch: epoch,
            mean_motion_dot: {
                let mut mstr: String = "0".to_owned();
                mstr.push_str(&line1[34..43]);
                let mut m: f64 = match mstr.parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse mean motion dot".to_string()),
                };
                if line1.chars().nth(33).unwrap() == '-' {
                    m = -1.0 * m;
                }
                m
            },
            mean_motion_dot_dot: {
                let mut mstr: String = "0.".to_owned();
                mstr.push_str(&line1[45..50]);
                mstr.push_str("E");
                mstr.push_str(&line1[50..53]);
                let mut m: f64 = match mstr.trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse mean motion dot dot".to_string()),
                };
                if line1.chars().nth(44).unwrap() == '-' {
                    m = -1.0 * m;
                }
                m
            },
            bstar: {
                let mut mstr: String = "0.".to_owned();
                mstr.push_str(&line1[54..59]);
                mstr.push_str("E");
                mstr.push_str(&line1[59..62]);
                let mut m: f64 = match mstr.trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse bstar (drag)".to_string()),
                };
                if line1.chars().nth(53).unwrap() == '-' {
                    m = -1.0 * m;
                }
                m
            },
            ephem_type: {
                match line1[62..63].trim().parse() {
                    Ok(y) => y,
                    Err(_) => 0,
                }
            },
            element_num: {
                match line1[64..68].trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse element number".to_string()),
                }
            },
            inclination: {
                match line2[8..16].trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse inclination".to_string()),
                }
            },
            raan: {
                match line2[17..25].trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse raan".to_string()),
                }
            },
            eccen: {
                let mut mstr: String = "0.".to_owned();
                mstr.push_str(&line2[26..33]);
                match mstr.trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse eccen".to_string()),
                }
            },
            arg_of_perigee: {
                match line2[34..42].trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse arg of perigee".to_string()),
                }
            },
            mean_anomaly: {
                match line2[42..51].trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse mean anomaly".to_string()),
                }
            },
            mean_motion: {
                match line2[52..63].trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse mean motion".to_string()),
                }
            },
            rev_num: {
                match line2[63..68].trim().parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse rev num".to_string()),
                }
            },
            satrec: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn testload() {
        let line1: &str = "1 26900U 01039A   06106.74503247  .00000045  00000-0  10000-3 0  8290";
        let line2: &str =
            "2 26900   0.0164 266.5378 0003319  86.1794 182.2590  1.00273847 16981   9300.";
        let line0: &str = "0 INTELSAT 902";

        //let line1: &str = "1 29238U 06022G   06177.28732010  .00766286  10823-4  13334-2 0   101";
        //let line2: &str = "2 29238  51.5595 213.7903 0202579  95.2503 267.9010 15.73823839  1061";
        //let line0: &str = "0 SL-12 DEB";

        match TLE::load_3line(&line0.to_string(), &line1.to_string(), &line2.to_string()) {
            Ok(t) => {
                assert_eq!(1, 1);
                println!("TLE = {:?}", t);
                println!("epoch = {}", t.epoch);
            }

            Err(s) => {
                println!("Err = \"{}\"", s);
                assert_eq!(1, 0);
            }
        }
    }
}
