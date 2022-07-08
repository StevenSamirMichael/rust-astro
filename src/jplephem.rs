//!
//! JPL Solar System Ephemerides
//!
//! # Introduction
//!
//! This module provides high-precision ephemerides
//! for bodies in the solar system, calculated from the
//! Jet Propulsion Laboratory (JPL) ephemeris data files
//!
//! ## Links
//!
//! Ephemerides filess can be found at:
//! <https://ssd.jpl.nasa.gov/ftp/eph/planets//>
//!
//!
//! ### Notes
//!
//! for little-endian systems, download from the "Linux" subdirectory
//! For big-endian systems, download from the "SunOS" subdirectory
//!

extern crate nalgebra;
use nalgebra::DMatrix;

use nalgebra as na;
pub type Vec3 = na::Vector3<f64>;
pub type Quat = na::UnitQuaternion<f64>;

use super::astrotime::{AstroTime, Scale};

/// Solar system bodies for which
/// JPL high-precision ephemeris can be queried
///
/// Coordinates origin is the solar system barycenter
///
/// ## Notes:
///  * Positions for all bodies are natively relative to solar system barycenter,
///    with exception of moon, which is computed in Geocentric system
///  * EMB (2) is the Earth-Moon barycenter
///  * The sun position is relative to the solar system barycenter
///    (it will be close to origin)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EphBody {
    // Mercury
    MERCURY = 0,
    // Venus
    VENUS = 1,
    // Earth-Moon Barycenter
    EMB = 2,
    // Mars
    MARS = 3,
    // Jupiter
    JUPITER = 4,
    // Saturn
    SATURN = 5,
    // Uranus
    URANUS = 6,
    // Neptune
    NEPTUNE = 7,
    // Pluto
    PLUTO = 8,
    // Moon (Geocentric)
    MOON = 9,
    // Sun
    SUN = 10,
}

impl TryFrom<i32> for EphBody {
    type Error = ();
    fn try_from(v: i32) -> Result<Self, Self::Error> {
        match v {
            x if x == EphBody::MERCURY as i32 => Ok(EphBody::MERCURY),
            x if x == EphBody::VENUS as i32 => Ok(EphBody::VENUS),
            x if x == EphBody::EMB as i32 => Ok(EphBody::EMB),
            x if x == EphBody::MARS as i32 => Ok(EphBody::MARS),
            x if x == EphBody::JUPITER as i32 => Ok(EphBody::JUPITER),
            x if x == EphBody::SATURN as i32 => Ok(EphBody::SATURN),
            x if x == EphBody::URANUS as i32 => Ok(EphBody::URANUS),
            x if x == EphBody::NEPTUNE as i32 => Ok(EphBody::NEPTUNE),
            x if x == EphBody::PLUTO as i32 => Ok(EphBody::PLUTO),
            x if x == EphBody::MOON as i32 => Ok(EphBody::MOON),
            x if x == EphBody::SUN as i32 => Ok(EphBody::SUN),
            _ => Err(()),
        }
    }
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

/// JPL Ephemeris Structure
///
/// included ephemerides and solar system constants loaded from the
/// JPL ephemerides file.
///
/// Also includes functions to compute heliocentric and geocentric
/// positions and velocities of solar system bodies as a function
/// of time
#[derive(Debug)]
pub struct JPLEphem {
    pub de_version: i32,
    /// Julian date of start of ephemerides database
    pub jd_start: f64,
    /// Julian date of end of ephemerides database
    pub jd_stop: f64,
    jd_step: f64,
    /// Length of 1 astronomical unit, km
    pub au: f64,
    /// Earth/Moon Ratio
    pub emrat: f64,

    // Offset lookup table
    ipt: [[usize; 3]; 15],
    consts: std::collections::HashMap<String, f64>,
    cheby: DMatrix<f64>,
}

impl JPLEphem {
    pub fn consts(&self, s: &String) -> Option<&f64> {
        self.consts.get(s)
    }

    /// Construct a JPL Ephemerides object from the provided binary data file
    ///
    ///
    /// # Return
    ///
    /// * Object holding all of the JPL ephemerides in memory that can be
    /// queried to find solar system body position in heliocentric or
    /// geocentric coordinate system as function of time
    ///
    /// # Example
    ///
    /// ```
    ///  // Note: filename can be full path. If not, file is
    ///  // assumed to be in the "datadir" directory
    ///  // (see submodule) that holds all the data files
    ///  // for the rest of the library
    ///
    /// use astro::jplephem::{JPLEphem, EphBody};
    /// use astro::astrotime::AstroTime;
    /// let jpl = JPLEphem::from_file("jplephem.440").unwrap();
    ///
    /// // Construct time: March 1 2021 12:00pm UTC
    /// let t = AstroTime::from_datetime(2021, 3, 1, 12, 0, 0.0);
    ///
    /// // Find geocentric moon position at this time
    /// let p = jpl.geocentric_body_pos(EphBody::MOON, &t).unwrap();
    /// ```
    ///
    pub fn from_file(fname: &str) -> Result<JPLEphem, Box<dyn std::error::Error>> {
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
                let mut v: DMatrix<f64> = DMatrix::repeat(ncoeff, nrecords, 0.0);

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

    // Optimized function for computing body position
    // (Matrix is allocated on stack, not heap)
    fn body_pos_optimized<const N: usize>(
        &self,
        body: EphBody,
        tm: &AstroTime,
    ) -> Result<Vec3, Box<dyn std::error::Error>> {
        // Terrestrial time
        let tt = tm.to_jd(Scale::TT);
        if (self.jd_start > tt) || (self.jd_stop < tt) {
            return Err(Box::new(InvalidTime));
        }

        // Get record index
        let t_int: f64 = (tt - self.jd_start) / self.jd_step;
        let int_num = t_int.floor() as i32;
        // Body index
        let bidx = body as usize;

        // # of coefficients and subintervals for this body
        let ncoeff = self.ipt[bidx][1];
        let nsubint = self.ipt[bidx][2];

        // Fractional way into step
        let t_int_2 = (t_int - int_num as f64) * nsubint as f64;
        let sub_int_num: usize = t_int_2.floor() as usize;
        // Scale from -1 to 1
        let t_seg = 2.0 * (t_int_2 - sub_int_num as f64) - 1.0;

        let offset0 = self.ipt[bidx][0] - 1 + sub_int_num * ncoeff * 3;

        let mut t = na::Vector::<f64, na::Const<N>, na::ArrayStorage<f64, N, 1>>::zeros();
        t[0] = 1.0;
        t[1] = t_seg;
        for j in 2..ncoeff {
            t[j] = 2.0 * t_seg * t[j - 1] - t[j - 2];
        }

        let mut pos: Vec3 = Vec3::zeros();
        for ix in 0..3 {
            let m = self
                .cheby
                .fixed_slice::<N, 1>(offset0 + N * ix, int_num as usize);
            pos[ix] = (m.transpose() * t)[(0, 0)];
        }

        Ok(pos * 1.0e3)
    }

    /// Return the position of the given body in the Heliocentric coordinate system
    ///
    /// # Inputs
    ///
    ///  * body - the solar system body for which to return position
    ///  * tm - The time at which to return position
    ///
    /// # Return
    ///
    ///    3-vector of cartesian Heliocentric position in meters
    ///
    ///
    /// # Notes:
    ///  * Positions for all bodies are natively relative to solar system barycenter,
    ///    with exception of moon, which is computed in Geocentric system
    ///  * EMB (2) is the Earth-Moon barycenter
    ///  * The sun position is relative to the solar system barycenter
    ///    (it will be close to origin)
    pub fn body_pos(
        &self,
        body: EphBody,
        tm: &AstroTime,
    ) -> Result<Vec3, Box<dyn std::error::Error>> {
        match self.ipt[body as usize][1] {
            6 => self.body_pos_optimized::<6>(body, tm),
            7 => self.body_pos_optimized::<7>(body, tm),
            8 => self.body_pos_optimized::<8>(body, tm),
            10 => self.body_pos_optimized::<10>(body, tm),
            11 => self.body_pos_optimized::<11>(body, tm),
            12 => self.body_pos_optimized::<12>(body, tm),
            13 => self.body_pos_optimized::<13>(body, tm),
            14 => self.body_pos_optimized::<14>(body, tm),
            _ => panic!("Invalid body"),
        }
    }
    /// Return the position & velocity the given body in the Heliocentric coordinate system
    ///
    /// # Inputs
    ///
    ///  * body - the solar system body for which to return position
    ///  * tm - The time at which to return position
    ///
    /// # Return
    ///
    ///  * Tuple with following values:
    ///    * 3-vector of cartesian Heliocentric position in meters
    ///    * 3-vector of cartesian Heliocentric velocity in meters / second
    ///
    ///
    /// # Notes:
    ///  * Positions for all bodies are natively relative to solar system barycenter,
    ///    with exception of moon, which is computed in Geocentric system
    ///  * EMB (2) is the Earth-Moon barycenter
    ///  * The sun position is relative to the solar system barycenter
    ///    (it will be close to origin)
    pub fn body_state(
        &self,
        body: EphBody,
        tm: &AstroTime,
    ) -> Result<(Vec3, Vec3), Box<dyn std::error::Error>> {
        match self.ipt[body as usize][1] {
            6 => self.body_state_optimized::<6>(body, tm),
            7 => self.body_state_optimized::<7>(body, tm),
            8 => self.body_state_optimized::<8>(body, tm),
            10 => self.body_state_optimized::<10>(body, tm),
            11 => self.body_state_optimized::<11>(body, tm),
            12 => self.body_state_optimized::<12>(body, tm),
            13 => self.body_state_optimized::<13>(body, tm),
            14 => self.body_state_optimized::<14>(body, tm),
            _ => panic!("Invalid body"),
        }
    }

    fn body_state_optimized<const N: usize>(
        &self,
        body: EphBody,
        tm: &AstroTime,
    ) -> Result<(Vec3, Vec3), Box<dyn std::error::Error>> {
        // Terrestrial time
        let tt = tm.to_jd(Scale::TT);
        if (self.jd_start > tt) || (self.jd_stop < tt) {
            return Err(Box::new(InvalidTime));
        }

        // Get record index
        let t_int: f64 = (tt - self.jd_start) / self.jd_step;
        let int_num = t_int.floor() as i32;
        // Body index
        let bidx = body as usize;

        // # of coefficients and subintervals for this body
        let ncoeff = self.ipt[bidx][1];
        let nsubint = self.ipt[bidx][2];

        // Fractional way into step
        let t_int_2 = (t_int - int_num as f64) * nsubint as f64;
        let sub_int_num: usize = t_int_2.floor() as usize;
        // Scale from -1 to 1
        let t_seg = 2.0 * (t_int_2 - sub_int_num as f64) - 1.0;

        let offset0 = self.ipt[bidx][0] - 1 + sub_int_num * ncoeff * 3;

        let mut t = na::Vector::<f64, na::Const<N>, na::ArrayStorage<f64, N, 1>>::zeros();
        let mut v = na::Vector::<f64, na::Const<N>, na::ArrayStorage<f64, N, 1>>::zeros();
        t[0] = 1.0;
        t[1] = t_seg;
        v[0] = 0.0;
        v[1] = 1.0;
        for j in 2..ncoeff {
            t[j] = 2.0 * t_seg * t[j - 1] - t[j - 2];
            v[j] = 2.0 * t_seg * v[j - 1] - v[j - 2] + 2.0 * t[j - 1];
        }

        let mut pos: Vec3 = Vec3::zeros();
        let mut vel: Vec3 = Vec3::zeros();
        for ix in 0..3 {
            let m = self
                .cheby
                .fixed_slice::<N, 1>(offset0 + N * ix, int_num as usize);
            pos[ix] = (m.transpose() * t)[(0, 0)];
            vel[ix] = (m.transpose() * v)[(0, 0)];
        }

        Ok((
            pos * 1.0e3,
            vel * 2.0e3 * nsubint as f64 / self.jd_step / 86400.0,
        ))
    }

    /// Return the position of the given body in the Geocentric coordinate system
    ///
    /// # Inputs
    ///
    ///  * body - the solar system body for which to return position
    ///  * tm - The time at which to return position
    ///
    /// # Return
    ///
    ///    3-vector of cartesian Geocentric position in meters
    ///
    pub fn geocentric_body_pos(
        &self,
        body: EphBody,
        tm: &AstroTime,
    ) -> Result<Vec3, Box<dyn std::error::Error>> {
        if body == EphBody::MOON {
            return self.body_pos(body, tm);
        } else {
            let emb: Vec3 = self.body_pos(EphBody::EMB, tm).unwrap();
            let moon: Vec3 = self.body_pos(EphBody::MOON, tm).unwrap();
            let b: Vec3 = self.body_pos(body, tm).unwrap();

            // Compute the position of the body relative to the Earth-moon
            // barycenter, then "correct" to Earth-center by accounting
            // for moon position and Earth/moon mass ratio
            Ok(b - emb + moon / (1.0 + self.emrat))
        }
    }

    /// Return the position and velocity of the given body in the Geocentric coordinate system
    ///
    /// # Inputs
    ///
    ///  * body - the solar system body for which to return position
    ///  * tm - The time at which to return position
    ///
    /// # Return
    ///
    ///   * Tuple with following elements:
    ///     * 3-vector of cartesian Geocentric position in meters
    ///     * 3-vector of cartesian Geocentric velocity in meters / second
    ///
    pub fn geocentric_body_state(
        &self,
        body: EphBody,
        tm: &AstroTime,
    ) -> Result<(Vec3, Vec3), Box<dyn std::error::Error>> {
        if body == EphBody::MOON {
            return self.body_state(body, tm);
        } else {
            let emb: Vec3 = self.body_pos(EphBody::EMB, tm).unwrap();
            let moon: Vec3 = self.body_pos(EphBody::MOON, tm).unwrap();
            let b: (Vec3, Vec3) = self.body_state(body, tm).unwrap();

            // Compute the position of the body relative to the Earth-moon
            // barycenter, then "correct" to Earth-center by accounting
            // for moon position and Earth/moon mass ratio
            Ok((b.0 - emb + moon / (1.0 + self.emrat), b.1))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Load the test vectors that come with the JPL ephemeris files
    /// and compare calculated positions to test vectors.
    #[test]
    fn testvecs() {
        let jpl = JPLEphem::from_file("jpleph.440").unwrap();

        // Read in test vectors from the NASA JPL web site
        // Note: this file will be stored in the "datadir" directory
        // permanently.
        use super::super::utils::dev::*;
        let lines =
            match lines_from_url("https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/testpo.440")
            {
                Ok(v) => (v),
                Err(e) => {
                    panic!("Cannot download test points: {}", e.to_string());
                }
            };
        for rline in lines.skip(14) {
            let line = match rline {
                Ok(v) => v,
                Err(_) => continue,
            };
            let s: Vec<&str> = line.split_whitespace().collect();
            assert!(s.len() >= 7);
            let jd: f64 = s[2].parse().unwrap();
            let tar: i32 = s[3].parse().unwrap();
            let src: i32 = s[4].parse().unwrap();
            let coord: usize = s[5].parse().unwrap();
            let truth: f64 = s[6].parse().unwrap();
            let tm = AstroTime::from_jd(jd, Scale::TT);
            if tar <= 10 && src <= 10 && coord <= 6 {
                let (mut tpos, mut tvel) = jpl
                    .geocentric_body_state(EphBody::try_from(tar - 1).unwrap(), &tm)
                    .unwrap();
                let (mut spos, mut svel) = jpl
                    .geocentric_body_state(EphBody::try_from(src - 1).unwrap(), &tm)
                    .unwrap();

                // in test vectors, index 3 is not EMB, but rather Earth
                // (this took me a long time to figure out...)
                if tar == 3 {
                    tpos = Vec3::zeros();
                    let (_mpos, mvel): (Vec3, Vec3) =
                        jpl.geocentric_body_state(EphBody::MOON, &tm).unwrap();
                    // Scale Earth velocity
                    tvel = tvel - mvel / (1.0 + jpl.emrat);
                }
                if src == 3 {
                    spos = Vec3::zeros();
                    let (_mpos, mvel): (Vec3, Vec3) =
                        jpl.geocentric_body_state(EphBody::MOON, &tm).unwrap();
                    //Scale Earth velocity
                    svel = svel - mvel / (1.0 + jpl.emrat);
                }
                if src == 10 {
                    // Compute moon velocity in heliocentric frame (not relative to Earth)
                    let (_embpos, embvel): (Vec3, Vec3) =
                        jpl.geocentric_body_state(EphBody::EMB, &tm).unwrap();
                    svel = svel + (embvel - svel / (1.0 + jpl.emrat));
                }
                if tar == 10 {
                    // Comput moon velocity in heliocentric frame (not relative to Earth)
                    let (_embpos, embvel): (Vec3, Vec3) =
                        jpl.geocentric_body_state(EphBody::EMB, &tm).unwrap();
                    tvel = tvel + (embvel - tvel / (1.0 + jpl.emrat));
                }

                // Comparing positions
                if coord <= 3 {
                    let calc = (tpos - spos)[coord - 1] / jpl.au / 1.0e3;
                    // These should be very exact
                    // Allow for errors of only ~ 1e-12
                    let maxerr = 1.0e-12;
                    let err = ((truth - calc) / truth).abs();
                    assert!(err < maxerr);
                }
                // Comparing velocities
                else {
                    let calc = (tvel - svel)[coord - 4] / jpl.au / 1.0e3 * 86400.0;
                    let maxerr: f64 = 1.0e-12;
                    let err: f64 = ((truth - calc) / truth).abs();
                    assert!(err < maxerr);
                }
            }
        }
    }

    #[test]
    fn load_test() {
        let jpl = JPLEphem::from_file("jpleph.440").unwrap();

        //let tm = &AstroTime::from_date(2010, 3, 1);
        let tm = AstroTime::from_jd(2451545.0, Scale::TT);
        println!("tm = {}", tm);
        println!("jd = {}", tm.to_jd(Scale::UTC));
        //let tm = &AstroTime::from_jd(2451545.0, Scale::UTC);
        let (pos, vel) = jpl.geocentric_body_state(EphBody::MOON, &tm).unwrap();
        println!("pos = {}", pos);
        println!("vel = {}", vel);
    }
}
