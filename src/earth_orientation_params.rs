use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;
use std::sync::RwLock;

use crate::utils::{astroerr, AstroResult};

use super::astrotime;
use crate::utils::datadir;
use crate::utils::download_file;
use lazy_static;
use std::sync::RwLock;

use crate::{astroerr, AstroResult};

#[derive(Debug)]
struct EOPEntry {
    mjd_utc: f64,
    xp: f64,
    yp: f64,
    dut1: f64,
    lod: f64,
}

fn load_eop_file(filename: Option<PathBuf>) -> Vec<EOPEntry> {
    let path: PathBuf = match filename {
        Some(pb) => pb,
        None => datadir::get()
            .unwrap_or(PathBuf::from("."))
            .join("finals2000A.all"),
    };

    if !path.is_file() {
        panic!(
            "Cannot open earth orientation parameters file: {}",
            path.to_str().unwrap()
        );
    }

    let file = match File::open(&path) {
        Err(why) => panic!("Couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };

    let mut eopvec = Vec::<EOPEntry>::new();
    for line in io::BufReader::new(file).lines() {
        match &line.unwrap() {
            v if v.len() < 100 => (),
            v if !v.is_ascii() => (),
            v if {
                let c: String = v.chars().skip(16).take(1).collect();
                c != "I" && c != "P"
            } =>
            {
                ()
            }
            v => {
                let mjd_str: String = v.chars().skip(7).take(8).collect();
                let xp_str: String = v.chars().skip(18).take(9).collect();
                let yp_str: String = v.chars().skip(37).take(9).collect();
                let dut1_str: String = v.chars().skip(58).take(10).collect();
                let lod_str: String = v.chars().skip(79).take(7).collect();

                eopvec.push(EOPEntry {
                    mjd_utc: {
                        match mjd_str.trim().parse() {
                            Ok(v) => v,
                            Err(_) => panic!("Could not extract mjd from file"),
                        }
                    },
                    xp: {
                        match xp_str.trim().parse() {
                            Ok(v) => v,
                            Err(_) => panic!("Could not extract x polar motion from file"),
                        }
                    },
                    yp: {
                        match yp_str.trim().parse() {
                            Ok(v) => v,
                            Err(_) => panic!("Could not extract y polar motion from file"),
                        }
                    },
                    dut1: {
                        match dut1_str.trim().parse() {
                            Ok(v) => v,
                            Err(_) => panic!("Could not extract dut1 from file"),
                        }
                    },
                    lod: lod_str.trim().parse().unwrap_or(0.0),
                })
            }
        }
    }
    eopvec
}

lazy_static::lazy_static! {
    static ref EOP_PARAMS: RwLock<Vec<EOPEntry>> = RwLock::new(load_eop_file(None));
}

/// Download new Earth Orientation Parameters file, and load it.
/// By default, tries to download into writeable directory where
/// data files for this crate are stored: crate::datadir::get()
pub fn update() -> AstroResult<()> {
    // Find writeabld data directory
    let d: Vec<PathBuf> = datadir::get_testdirs()
        .into_iter()
        .filter(|x| x.is_dir())
        .filter(|x| x.metadata().unwrap().permissions().readonly() == false)
        .collect();
    if d.len() == 0 {
        return astroerr!("Cannot find writable data directory");
    }

    // Download most-recent EOP
    let url = "https://datacenter.iers.org/data/9/finals2000A.all";
    download_file(url, &d[0], true)?;

    // Re-load the params
    *EOP_PARAMS.write().unwrap() = load_eop_file(None);

    Ok(())
}

lazy_static::lazy_static! {
    static ref EOP_PARAMS : RwLock<EOPData> = RwLock::new(load_eop_data(None));
}

///
/// Download updated Earth Orientation Parameters
/// to provided dita directory.
///
/// These should be updated daily
///
///
pub fn update() -> AstroResult<()> {
    let path = datadir::get()?;
    if path.metadata()?.permissions().readonly() {
        let s = path.to_str().unwrap();
        return astroerr!("Cannot write to data directory: {}", s);
    }
    let url = "https://celestrak.org/SpaceData/sw19571001.txt";
    download_file(url, &path, true)?;

    let mut eopguard = EOP_PARAMS.write().unwrap();
    *eopguard = load_eop_data(None);

    Ok(())
}

///
/// Get Earth orientation parameters given UTC modified Julian date
///
pub fn get_from_mjd_utc(mjd_utc: f64) -> Option<[f64; 4]> {
    let eop = EOP_PARAMS.read().unwrap();

    let idx = eop.iter().position(|x| x.mjd_utc > mjd_utc);
    match idx {
        None => None,
        Some(v) => {
            if v == 0 as usize {
                return None;
            }
            // Linear interpolation
            let g1: f64 = (mjd_utc - eop[v - 1].mjd_utc) / (eop[v].mjd_utc - eop[v - 1].mjd_utc);
            let g0: f64 = 1.0 - g1;
            let v1: &EOPEntry = &eop[v];
            let v0: &EOPEntry = &eop[v - 1];
            Some([
                g0 * v0.dut1 + g1 * v1.dut1,
                g0 * v0.xp + g1 * v1.xp,
                g0 * v0.yp + g1 * v1.yp,
                g0 * v0.lod + g1 * v1.lod,
            ])
        }
    }
}

///
/// Get Earth Orientation Parameters at given instant
///
/// # Arguments:
///
/// * tm: Instant at which to query parameters
///
/// # Return:
///
/// * Vector [f64; 4] with following elements:
///   * 0 : (UT1 - UTC) in seconds
///   * 1 : X polar motion in arcsecs
///   * 2 : Y polar motion in arcsecs
///   * 3 : LOD: instantaneous rate of change in (UT1-UTC), msec/day
///
#[inline]
pub fn get(tm: &astrotime::AstroTime) -> Option<[f64; 4]> {
    get_from_mjd_utc(tm.to_mjd(astrotime::Scale::UTC))
}

#[cfg(test)]
mod tests {

    use super::*;

    /// Check that data is loaded
    #[test]
    fn loaded() {
        assert_eq!(EOP_PARAMS.read().unwrap()[0].mjd_utc >= 0.0, true);
    }

    /// Check value against manual value from file
    #[test]
    fn checkval() {
        let v = get_from_mjd_utc(59464.00).unwrap();
        const TRUTH: [f64; 4] = [-0.1145681, 0.241135, 0.317269, -0.2418];
        for it in v.iter().zip(TRUTH.iter()) {
            let (a, b) = it;
            assert!(((a - b) / b).abs() < 1.0e-5);
        }
    }

    /// Check interpolation between points
    #[test]
    fn checkinterp() {
        let mjd0: f64 = 57909.00;
        const TRUTH0: [f64; 4] = [0.3754796, 0.10268, 0.458479, 1.1596];
        const TRUTH1: [f64; 4] = [0.3743671, 0.10402, 0.458388, 1.0584];
        for x in 0..101 {
            let dt: f64 = x as f64 / 100.0;
            let vt = get_from_mjd_utc(mjd0 + dt).unwrap();
            let g0: f64 = 1.0 - dt;
            let g1: f64 = dt;
            for it in vt.iter().zip(TRUTH0.iter().zip(TRUTH1.iter())) {
                let (v, (v0, v1)) = it;
                let vtest: f64 = g0 * v0 + g1 * v1;
                assert!(((v - vtest) / v).abs() < 1.0e-5);
            }
        }
    }
}
