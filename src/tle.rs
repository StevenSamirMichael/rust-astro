use super::astrotime::AstroTime;

#[derive(PartialEq, PartialOrd, Clone)]
pub struct TLE {
    pub name: String,
    pub intl_desig: String,
    pub sat_num: i32,
    pub desig_year: i32,
    pub desig_launch: i32,
    pub desig_piece: String,
    pub epoch: AstroTime,
    pub mean_motion_dot: f64,
    pub mean_motion_dot_dot: f64,
    pub bstar: f64,
    pub ephem_type: u8,
    pub element_num: i32,
    pub inclination: f64,
    pub raan: f64,
    pub eccen: f64,
    pub arg_of_perigee: f64,
    pub mean_anomaly: f64,
    pub mean_motion: f64,
    pub rev_num: i32,

    need_reinit: bool,
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
            need_reinit: true,
        }
    }

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
        let mut year: u32 = match line1[18..21].parse() {
            Ok(y) => y,
            Err(_) => return Err("Could not parse year".to_string()),
        };
        if year > 57 {
            year += 1900;
        } else {
            year += 2000;
        }
        let day_of_year: f64 = match line1[20..33].parse() {
            Ok(y) => y,
            Err(_) => return Err("Could not parse day of year".to_string()),
        };
        let i_day_of_year: u32 = day_of_year.floor() as u32;
        let frac_of_day: f64 = day_of_year - (i_day_of_year as f64);
        let epoch: AstroTime =
            AstroTime::from_datetime(year, 1, i_day_of_year, 0, 0, frac_of_day * 86400.0);

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
                match line1[9..12].trim().parse() {
                    Ok(l) => l,
                    Err(_) => return Err("Could not parse desig year".to_string()),
                }
            },
            desig_launch: {
                match line1[11..14].trim().parse() {
                    Ok(l) => l,
                    Err(_) => return Err("Could not parse desig_launch".to_string()),
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
                let mut m: f64 = match mstr.parse() {
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
                let mut m: f64 = match mstr.parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse bstar (drag)".to_string()),
                };
                if line1.chars().nth(53).unwrap() == '-' {
                    m = -1.0 * m;
                }
                m
            },
            ephem_type: {
                match line1[62..63].parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse ephem type".to_string()),
                }
            },
            element_num: {
                match line1[64..68].parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse element number".to_string()),
                }
            },
            inclination: {
                match line2[8..16].parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse inclination".to_string()),
                }
            },
            raan: {
                match line2[17..25].parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse raan".to_string()),
                }
            },
            eccen: {
                let mut mstr: String = "0.".to_owned();
                mstr.push_str(&line2[26..33]);
                match mstr.parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse eccen".to_string()),
                }
            },
            arg_of_perigee: {
                match line2[34..42].parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse arg of perigee".to_string()),
                }
            },
            mean_anomaly: {
                match line2[42..51].parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse mean anomaly".to_string()),
                }
            },
            mean_motion: {
                match line2[52..63].parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse mean motion".to_string()),
                }
            },
            rev_num: {
                match line2[63..68].parse() {
                    Ok(y) => y,
                    Err(_) => return Err("Could not parse rev num".to_string()),
                }
            },
            need_reinit: true,
        })
    }
}
