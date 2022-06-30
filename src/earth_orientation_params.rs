use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;

use super::astrotime;
use super::datadir;
use lazy_static;

struct EOPEntry {
    mjd_utc: f64,
    xp: f64,
    yp: f64,
    dut1: f64,
    lod: f64,
}

lazy_static::lazy_static! {
    static ref EOP_PARAMS: Vec<EOPEntry> = {
        let path = datadir::get()
            .unwrap_or(PathBuf::from("."))
            .join("finals2000A.all");
        if !path.is_file() {
            panic!("Cannot open earth orientation parameters file: {}", path.to_str().unwrap());
        }

        let file = match File::open(&path) {
            Err(why) => panic!("Couldn't open {}: {}", path.display(), why),
            Ok(file) => file,
        };

        let mut eopvec = Vec::<EOPEntry>::new();
        for line in io::BufReader::new(file).lines() {
            match &line {
                Ok(s) => match &s {
                    v if v.len() < 100 => (),
                    v if !v.is_ascii() => (),
                    v => {
                        let mjd_str: String =
                            v.chars().skip(7).take(8).collect();
                        let xp_str: String =
                            v.chars().skip(18).take(9).collect();

                        let yp_str: String =
                            v.chars().skip(37).take(9).collect();

                        let dut1_str: String =
                            v.chars().skip(58).take(10).collect();

                        let lod_str: String =
                            v.chars().skip(79).take(7).collect();

                        eopvec.push(EOPEntry {
                            mjd_utc: {
                                match mjd_str.parse() {
                                Ok(v)=> v,
                                Err(_) => panic!("Could not extract mjd from file")
                               }
                            },
                            xp: {
                                match xp_str.parse() {
                                Ok(v)=>v,
                                Err(_)=>panic!("Could not extract x polar motion from file")
                                }
                            },
                            yp: {
                                match yp_str.parse() {
                                Ok(v)=>v,
                                Err(_)=>panic!("Could not extract y polar motion from file")
                                }
                            },
                            dut1: {
                                match dut1_str.parse() {
                                    Ok(v)=>v,
                                    Err(_)=>panic!("Could not extract dut1 from file")
                                }
                            },
                            lod: lod_str.parse().unwrap_or(0.0)
                        })
                    }
                },
                Err(_) => (),
            }
        }
        eopvec
    };
}

pub fn get_from_mjd_utc(mjd_utc: f64) -> Option<(f64, f64, f64, f64)> {
    let idx = EOP_PARAMS.iter().position(|x| x.mjd_utc > mjd_utc);
    match idx {
        None => None,
        Some(v) => {
            if v == 0 as usize {
                return None;
            }
            // Linear interpolation
            let g1: f64 = (mjd_utc - EOP_PARAMS[v - 1].mjd_utc)
                / (EOP_PARAMS[v].mjd_utc - EOP_PARAMS[v - 1].mjd_utc);
            let g0: f64 = 1.0 - g1;
            let v1: &EOPEntry = &EOP_PARAMS[v];
            let v0: &EOPEntry = &EOP_PARAMS[v - 1];
            Some((
                g0 * v0.dut1 + g1 * v1.dut1,
                g0 * v0.xp + g1 * v1.xp,
                g0 * v0.yp + g1 * v1.yp,
                g0 * v0.lod + g1 * v1.lod,
            ))
        }
    }
}

#[inline]
pub fn get(tm: &astrotime::AstroTime) -> Option<(f64, f64, f64, f64)> {
    get_from_mjd_utc(tm.to_mjd(astrotime::Scale::UTC))
}
