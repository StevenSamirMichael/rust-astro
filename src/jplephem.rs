use nalgebra::DMatrix;

use nalgebra as na;
pub type Vec3 = na::Vector3<f64>;
pub type Quat = na::UnitQuaternion<f64>;

use super::astrotime::{AstroTime, Scale};

/// Solar system bodies for which
/// JPL ephemeris can be queried
///
/// Coordinates origin is the solar system barycenter
///
/// EMB (2) is the Earth-Moon barycenter
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EphBody {
    MERCURY = 0,
    VENUS = 1,
    EMB = 2,
    MARS = 3,
    JUPITER = 4,
    SATURN = 5,
    URANUS = 6,
    NEPTUNE = 7,
    PLUTO = 8,
    MOON = 9,
    SUN = 10,
}

#[derive(Debug, Clone)]
struct InvalidSize;
impl std::error::Error for InvalidSize {}
impl std::fmt::Display for InvalidSize {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Invalid file size")
    }
}

#[derive(Debug, Clone)]
struct InvalidTime;
impl std::error::Error for InvalidTime {}
impl std::fmt::Display for InvalidTime {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Invalid Time for Body Query")
    }
}

#[derive(Debug)]
struct JPLEphem {
    de_version: i32,
    jd_start: f64,
    jd_stop: f64,
    jd_step: f64,
    au: f64,
    emrat: f64,
    ipt: [[usize; 3]; 15],
    consts: std::collections::HashMap<String, f64>,
    cheby: DMatrix<f64>,
}

pub fn emrat() -> f64 {
    JPLEPHEM.emrat
}

pub fn au() -> f64 {
    JPLEPHEM.au
}

pub fn consts(s: &String) -> Option<&f64> {
    JPLEPHEM.consts.get(s)
}

pub fn de_version() -> i32 {
    JPLEPHEM.de_version
}

fn load_ephemeris_file(fname: &str) -> Result<JPLEphem, Box<dyn std::error::Error>> {
    use super::datadir;
    use std::collections::HashMap;
    use std::path::PathBuf;

    // Dimensions of ephemeris for given index
    fn dimension(idx: usize) -> usize {
        if idx == 11 {
            2
        } else if idx == 14 {
            1
        } else {
            3
        }
    }

    // Open the file
    let path = datadir::get().unwrap_or(PathBuf::from(".")).join(fname);
    if !path.is_file() {
        panic!("Cannot open JPL Ephemeris file");
    }

    // Read in bytes
    let raw = std::fs::read(path)?;
    let title: &str = &std::str::from_utf8(&raw[0..84])?;

    // Get version
    let de_version: i32 = title[26..29].parse()?;

    let jd_start = f64::from_le_bytes(raw[2652..2660].try_into()?);
    let jd_stop: f64 = f64::from_le_bytes(raw[2660..2668].try_into()?);
    let jd_step: f64 = f64::from_le_bytes(raw[2668..2676].try_into()?);
    let n_con: i32 = i32::from_le_bytes(raw[2676..2680].try_into()?);
    let au: f64 = f64::from_le_bytes(raw[2680..2688].try_into()?);
    let emrat: f64 = f64::from_le_bytes(raw[2688..2696].try_into()?);

    // Get table
    let ipt: [[usize; 3]; 15] = {
        let mut ipt: [[usize; 3]; 15] = [[0, 0, 0]; 15];
        let mut idx = 2696;
        for ix in 0..15 {
            for iy in 0..3 {
                ipt[ix][iy] = u32::from_le_bytes(raw[idx..(idx + 4)].try_into()?) as usize;
                idx = idx + 4;
            }
        }

        ipt[12][0] = ipt[12][1];
        ipt[12][1] = ipt[12][2];
        ipt[12][2] = ipt[13][0];

        if de_version > 430 && n_con != 400 {
            if n_con > 400 {
                let idx = ((n_con - 400) * 6) as usize;
                ipt[13][0] = u32::from_le_bytes(raw[idx..(idx + 4)].try_into()?) as usize;
            } else {
                ipt[13][0] = 1 as usize;
            }
        }

        // Check for garbage data not populated in earlier files
        if ipt[13][0] != (ipt[12][0] + ipt[12][1] * ipt[12][2] * 3)
            || ipt[14][0] != (ipt[13][0] + ipt[13][1] * ipt[13][2] * 3)
        {
            for ix in 13..15 {
                for iy in 0..3 {
                    ipt[ix][iy] = 0;
                }
            }
        }
        ipt
    };

    // Kernel size
    let kernel_size: usize = {
        let mut ks: usize = 4;
        for ix in 0..15 {
            ks = ks + 2 * ipt[ix][1] * ipt[ix][2] * dimension(ix)
        }
        ks
    };

    Ok(JPLEphem {
        de_version: de_version,
        jd_start: jd_start,
        jd_stop: jd_stop,
        jd_step: jd_step,
        au: au,
        emrat: emrat,
        ipt: ipt,
        consts: {
            let mut hm = HashMap::new();

            // Read in constants defined in file
            for ix in 0..n_con {
                let sidx: usize = kernel_size * 4 + ix as usize * 8;
                let eidx: usize = sidx + 8;
                let val: f64 = f64::from_le_bytes(raw[sidx..eidx].try_into()?);

                let mut stridx: usize = (84 * 3 + ix * 6) as usize;
                // different loc if constants >= 400
                if ix >= 400 {
                    stridx = (84 * 3 + 400 * 6 + 5 * 8 + 41 * 4 + ix * 6) as usize;
                }
                let s = String::from_utf8(raw[stridx..(stridx + 6)].to_vec())?;

                hm.insert(String::from(s.trim()), val);
            }
            hm
        },
        cheby: {
            // we are going to do this unsafe since I can't find a
            // fast way to do it otherwise
            let ncoeff: usize = (kernel_size / 2) as usize;
            let nrecords = ((jd_stop - jd_start) / jd_step) as usize;
            let record_size = (kernel_size * 4) as usize;
            let mut v: DMatrix<f64> = DMatrix::from_element(ncoeff, nrecords, 0.0);
            if raw.len() < record_size * 2 + ncoeff * nrecords * 8 {
                return Err(Box::new(InvalidSize));
            }

            unsafe {
                std::ptr::copy_nonoverlapping(
                    raw.as_ptr().offset((record_size * 2) as isize) as *const f64,
                    v.as_mut_ptr() as *mut f64,
                    ncoeff * nrecords,
                );
            }
            v
        },
    })
}

lazy_static::lazy_static! {
    #[derive(Debug)]
    static ref JPLEPHEM: JPLEphem =
        load_ephemeris_file(&"jpleph.440").unwrap_or_else(|e| {panic!("Error: {}", e);});
}

pub fn body_pos_optimized<const N: usize>(
    body: EphBody,
    tm: &AstroTime,
) -> Result<Vec3, Box<dyn std::error::Error>> {
    // Terrestrial time
    let tt = tm.to_jd(Scale::TT);
    if (JPLEPHEM.jd_start > tt) || (JPLEPHEM.jd_stop < tt) {
        return Err(Box::new(InvalidTime));
    }

    // Get record intex
    let t_int: f64 = (tt - JPLEPHEM.jd_start) / JPLEPHEM.jd_step;
    let int_num = t_int.floor() as i32;
    // Body index
    let bidx = body as usize;

    // # of coefficients and subintervals for this body
    let ncoeff = JPLEPHEM.ipt[bidx][1];
    let nsubint = JPLEPHEM.ipt[bidx][2];

    // Fractional way into step
    let t_int_2 = (t_int - int_num as f64) * nsubint as f64;
    let sub_int_num: usize = t_int_2.floor() as usize;
    // Scale from -1 to 1
    let t_seg = 2.0 * (t_int_2 - sub_int_num as f64) - 1.0;

    let offset0 = JPLEPHEM.ipt[bidx][0] - 1 + sub_int_num * ncoeff * 3;

    let mut t = na::Vector::<f64, na::Const<N>, na::ArrayStorage<f64, N, 1>>::zeros();
    t[0] = 1.0;
    t[1] = t_seg;
    for j in 2..ncoeff {
        t[j] = 2.0 * t_seg * t[j - 1] - t[j - 2];
    }

    let mut pos: Vec3 = Vec3::zeros();
    for ix in 0..3 {
        let m = JPLEPHEM
            .cheby
            .fixed_slice::<N, 1>(offset0 + N * ix, int_num as usize);
        pos[ix] = (m.transpose() * t)[(0, 0)];
    }

    Ok(pos)
}

pub fn body_pos(body: EphBody, tm: &AstroTime) -> Result<Vec3, Box<dyn std::error::Error>> {
    match JPLEPHEM.ipt[body as usize][1] {
        6 => body_pos_optimized::<6>(body, tm),
        7 => body_pos_optimized::<7>(body, tm),
        8 => body_pos_optimized::<8>(body, tm),
        10 => body_pos_optimized::<10>(body, tm),
        11 => body_pos_optimized::<11>(body, tm),
        12 => body_pos_optimized::<12>(body, tm),
        13 => body_pos_optimized::<13>(body, tm),
        _ => panic!("Invalid body"),
    }
}

pub fn geocentric_body_pos(
    body: EphBody,
    tm: &AstroTime,
) -> Result<Vec3, Box<dyn std::error::Error>> {
    if body == EphBody::MOON {
        return body_pos(body, tm);
    } else {
        let emb: Vec3 = body_pos(EphBody::EMB, tm).unwrap();
        let moon: Vec3 = body_pos(EphBody::MOON, tm).unwrap();
        let b: Vec3 = body_pos(body, tm).unwrap();

        Ok(b - emb + JPLEPHEM.emrat * moon)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_test() {
        //let s = load_ephemeris_file("jpleph.440");
        //println!("s = {:?}", s.unwrap().cheby);
        let s = body_pos(EphBody::MOON, &AstroTime::from_jd(2451545.0, Scale::UTC)).unwrap();
        println!("s = {} {} {}", s[0], s[1], s[2]);
        println!("snorm = {} km", s.norm());

        assert!(1 == 1);
    }
}
