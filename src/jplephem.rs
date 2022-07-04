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
    n_con: i32,
    au: f64,
    emrat: f64,
    kernel_size: u32,
    record_size: u32,
    ncoeff: u32,
    ipt: [[u32; 3]; 15],
    consts: std::collections::HashMap<String, f64>,
    cheby: DMatrix<f64>,
}

fn load_ephemeris_file(fname: &str) -> Result<JPLEphem, Box<dyn std::error::Error>> {
    use super::datadir;
    use std::collections::HashMap;
    use std::path::PathBuf;

    // Dimensions of ephemeris for given index
    fn dimension(idx: usize) -> u32 {
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

    // Get number of constants
    let n_con: i32 = i32::from_le_bytes(raw[2676..2680].try_into()?);

    // Get table
    let ipt: [[u32; 3]; 15] = {
        let mut ipt: [[u32; 3]; 15] = [[0, 0, 0]; 15];
        let mut idx = 2696;
        for ix in 0..15 {
            for iy in 0..3 {
                ipt[ix][iy] = u32::from_le_bytes(raw[idx..(idx + 4)].try_into()?);
                idx = idx + 4;
            }
        }

        ipt[12][0] = ipt[12][1];
        ipt[12][1] = ipt[12][2];
        ipt[12][2] = ipt[13][0];

        if de_version > 430 && n_con != 400 {
            if n_con > 400 {
                let idx = ((n_con - 400) * 6) as usize;
                ipt[13][0] = u32::from_le_bytes(raw[idx..(idx + 4)].try_into()?);
            } else {
                ipt[13][0] = 1 as u32;
            }
        }
        // Check for garbage data not populated in earlier files
        if ipt[13][0] != (ipt[12][0] + ipt[12][1] + ipt[12][2]) * 3
            || ipt[14][0] != (ipt[13][0] + ipt[13][1] + ipt[13][2]) * 3
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
    let kernel_size: u32 = {
        let mut ks: u32 = 4;
        for ix in 0..15 {
            ks = ks + 2 * ipt[ix][1] * ipt[ix][2] * dimension(ix)
        }
        ks
    };

    // Julian date start, stop, & step size
    let jd_start = f64::from_le_bytes(raw[2652..2660].try_into()?);
    let jd_stop: f64 = f64::from_le_bytes(raw[2660..2668].try_into()?);
    let jd_step: f64 = f64::from_le_bytes(raw[2668..2676].try_into()?);

    Ok(JPLEphem {
        de_version: de_version,
        jd_start: jd_start,
        jd_stop: jd_stop,
        jd_step: jd_step,
        n_con: n_con,
        au: f64::from_le_bytes(raw[2680..2688].try_into()?),
        emrat: f64::from_le_bytes(raw[2688..2696].try_into()?),
        ipt: ipt,
        kernel_size: kernel_size,
        record_size: kernel_size * 4,
        ncoeff: kernel_size / 2,
        consts: {
            let mut hm = HashMap::new();

            // Read in constants defined in file
            for ix in 0..n_con {
                let sidx: usize = (kernel_size * 4 + ix as u32 * 8) as usize;
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
                std::ptr::copy(
                    raw.as_ptr().offset((record_size * 2) as isize) as *mut f64,
                    v.as_mut_ptr(),
                    ncoeff * nrecords,
                )
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

pub fn body_pos<const N: usize>(
    body: EphBody,
    tm: AstroTime,
) -> Result<Vec3, Box<dyn std::error::Error>> {
    let tt = tm.to_jd(Scale::TT);
    if (JPLEPHEM.jd_start > tt) || (JPLEPHEM.jd_stop < tt) {
        return Err(Box::new(InvalidTime));
    }

    let t_int: f64 = (tt - JPLEPHEM.jd_start) / JPLEPHEM.jd_step;
    let int_num = t_int.floor() as i32;
    let bidx = body as usize;

    let ncoeff = JPLEPHEM.ipt[bidx][1];
    let nsubint = JPLEPHEM.ipt[bidx][2];

    let t_int_2 = (t_int - int_num as f64) * nsubint as f64;
    let sub_int_num: u32 = t_int_2.floor() as u32;
    let t_seg = 2.0 * (t_int_2 - sub_int_num as f64) - 1.0;

    let offset0 = JPLEPHEM.ipt[bidx][0] - 1 + sub_int_num * ncoeff * 3;

    let mut t = na::Vector::<f64, na::Const<N>, na::ArrayStorage<f64, N, 1>>::zeros();
    t[0] = 1.0;
    t[1] = t_seg;
    println!("ncoeff = {}", ncoeff);
    for j in 2..ncoeff as usize {
        t[j] = 2.0 * t_seg + t[j - 1] - t[j - 2];
    }

    println!("t_sg = {}", t_seg);
    println!("offset0 = {}", offset0);
    println!("int_num = {}", int_num);

    let v: Vec3 = JPLEPHEM
        .cheby
        .fixed_columns::<3>(offset0 as usize)
        .fixed_rows::<N>(int_num as usize)
        .transpose()
        * t;

    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_test() {
        //let s = load_ephemeris_file("jpleph.440");
        //println!("s = {:?}", s.unwrap().cheby);
        println!(
            "s = {:?}",
            body_pos::<13>(EphBody::EMB, AstroTime::from_jd(2451545.0, Scale::TT))
        );
        assert!(1 == 1);
    }
}
