use super::sgp4_lowlevel::sgp4_lowlevel; // propagator
use super::sgp4init::sgp4init;

use crate::astrotime::{AstroTime, Scale};
use crate::tle::TLE;
use nalgebra::{Const, Dyn, OMatrix};

type StateArr = OMatrix<f64, Const<3>, Dyn>;
pub type SGP4State = (StateArr, StateArr);
pub type SGP4Result = Result<SGP4State, (i32, String)>;

use std::f64::consts::PI;

use super::{GravConst, OpsMode};

#[inline]
pub fn sgp4(tle: &mut TLE, tm: &[AstroTime]) -> SGP4Result {
    sgp4_full(tle, tm, GravConst::WGS84, OpsMode::IMPROVED)
}

const SGP4_ERRS: [&str; 7] = [
    "Success",
    "Mean elements, ecc > 1.0 or ecc < -0.001 or a  < 0.95 Earth radius",
    "Mean motion < 0.0",
    "Pert. Elements: ecc < 0.0 or ecc > 1.0",
    "Semi-latus rectum < 0.0",
    "Epoch elements are sub-orbital",
    "Satellite has decayed",
];

///
/// Run Simplified General Perturbations (SGP)-4 propagator on
/// Two-Line Element Set to
/// output satellite position and velocity at given time
/// in the "TEME" coordinate system
///
/// A detailed description is
/// [here](https://celestrak.org/publications/AIAA/2008-6770/AIAA-2008-6770.pdf)
///
///
/// # Arguments
///
/// * tle: The TLE on which top operate.  Note that a mutable reference
///   is passed, as SGP-4 metadata is stored after each propagation
/// * tm: The time at which to compute position and velocity
///
/// # Return
///
/// Result object containing either an OK value containing a tuple with
/// position (m) and velocity (m/s) 3-vectors or an Err value containing
/// a tuple with error code and error string
///
/// # Example
///
/// ```
/// // Compute the Geodetic position of a satellite at
/// // the TLE epoch time
///     
/// use astro::TLE;
/// use astro::sgp4::{sgp4_full, GravConst, OpsMode};
/// use astro::coordconversion::qteme2itrf;
/// use astro::itrfcoord::ITRFCoord;
///
/// let line0: &str = "0 INTELSAT 902";
/// let line1: &str = "1 26900U 01039A   06106.74503247  .00000045  00000-0  10000-3 0  8290";
/// let line2: &str = "2 26900   0.0164 266.5378 0003319  86.1794 182.2590  1.00273847 16981   9300.";
/// let mut tle = TLE::load_3line(&line0.to_string(),
///     &line1.to_string(),
///     &line2.to_string()
///     ).unwrap();
///
/// let tm = tle.epoch;
/// let (pteme, vteme) = sgp4_full(&mut tle,
///     &tm,
///     GravConst::WGS84,
///     OpsMode::IMPROVED
///     ).unwrap();
///
/// let pitrf = qteme2itrf(&tm) * pteme;
/// let itrf = ITRFCoord::from_vec(pitrf.into());
///
/// ```
///
pub fn sgp4_full<'a>(
    tle: &mut TLE,
    tm: &[AstroTime],
    gravconst: GravConst,
    opsmode: OpsMode,
) -> SGP4Result {
    const TWOPI: f64 = PI * 2.0;

    if tle.satrec.is_none() {
        let no = tle.mean_motion / (1440.0 / TWOPI);
        let bstar = tle.bstar;
        let ndot = tle.mean_motion_dot / (1440.0 * 1440.0 / TWOPI);
        let nddot = tle.mean_motion_dot_dot / (1440.0 * 1440.0 * 1440.0 / TWOPI);
        let inclo = tle.inclination * PI / 180.0;
        let nodeo = tle.raan * PI / 180.0;
        let argpo = tle.arg_of_perigee * PI / 180.0;
        let mo = tle.mean_anomaly * PI / 180.0;
        let ecco = tle.eccen;
        let jdsatepoch = tle.epoch.to_jd(Scale::UTC);

        match sgp4init(
            gravconst,
            opsmode,
            &"satno",
            jdsatepoch - 2433281.5,
            bstar,
            ndot,
            nddot,
            ecco,
            argpo,
            inclo,
            mo,
            no,
            nodeo,
        ) {
            Ok(sr) => tle.satrec = Some(sr),
            Err(e) => return Err((e, String::from(SGP4_ERRS[e as usize]))),
        }
    }

    let mut s = tle.satrec.as_mut().unwrap();

    let n = tm.len();
    let mut rarr = StateArr::zeros(n);
    let mut varr = StateArr::zeros(n);
    for (pos, thetime) in tm.iter().enumerate() {
        let tsince = (*thetime - tle.epoch) * 1440.0;

        match sgp4_lowlevel(&mut s, tsince) {
            Ok((r, v)) => {
                rarr.index_mut((.., pos)).copy_from_slice(&r);
                varr.index_mut((.., pos)).copy_from_slice(&v);
            }
            Err(e) => return Err((e, String::from(SGP4_ERRS[e as usize]))),
        }
    }
    Ok((rarr * 1.0e3, varr * 1.0e3))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tle::TLE;
    use crate::utils::test;
    use std::io::BufRead;

    #[test]
    fn testsgp4() {
        let line1: &str = "1 26900U 01039A   06106.74503247  .00000045  00000-0  10000-3 0  8290";
        let line2: &str =
            "2 26900   0.0164 266.5378 0003319  86.1794 182.2590  1.00273847 16981   9300.";
        let line0: &str = "0 INTELSAT 902";

        //let line1: &str = "1 29238U 06022G   06177.28732010  .00766286  10823-4  13334-2 0   101";
        //let line2: &str = "2 29238  51.5595 213.7903 0202579  95.2503 267.9010 15.73823839  1061";
        //let line0: &str = "0 SL-12 DEB";
        let mut tle =
            TLE::load_3line(&line0.to_string(), &line1.to_string(), &line2.to_string()).unwrap();
        let tm = tle.epoch;

        match sgp4(&mut tle, &[tm]) {
            Ok((pos, vel)) => {
                println!("pos = {}", pos);
                println!("vel = {}", vel);
            }
            Err(e) => {
                panic!("Error running sgp4: \"{}\"", e.1);
            }
        }
        println!("module path = {}", std::module_path!());
    }

    #[test]
    fn vallado_testvecs() {
        let testdir = test::get_project_root()
            .unwrap()
            .join("testdata")
            .join("vallado")
            .join("TestSGP4")
            .join("TestSGP4");
        if !testdir.is_dir() {
            panic!(
                "Required test directory \"{}\" does not exist
                try running \"getfiles.sh\" script in $(crateroot)/testdata",
                testdir.to_string_lossy()
            );
        }
        let tlefile = testdir.join("SGP4-VER.TLE");
        let f = match std::fs::File::open(&tlefile) {
            Err(why) => panic!("Could not open {}: {}", tlefile.display(), why),
            Ok(file) => file,
        };
        let buf = std::io::BufReader::new(f);
        let lines: Vec<String> = buf.lines().map(|l| l.unwrap()).collect();
        let tles = TLE::from_lines(&lines).unwrap();
        assert!(tles.len() > 5);

        for mut tle in tles {
            let fname = format!("{:05}.e", tle.sat_num);

            let fh = testdir.join(fname);
            let ftle = match std::fs::File::open(&fh) {
                Err(why) => panic!("Could not open {}: {}", fh.display(), why),
                Ok(file) => file,
            };
            for line in std::io::BufReader::new(ftle).lines() {
                let maxposerr = 1.0e-5;
                let mut maxvelerr = 1.0e-5;

                let testvec: Vec<f64> = line
                    .unwrap()
                    .trim()
                    .split_whitespace()
                    .map(|x| match x.parse() {
                        Ok(v) => v,
                        Err(_) => -1.0,
                    })
                    .collect();
                if testvec.len() < 7 {
                    continue;
                }
                if testvec[0] < 0.0 {
                    continue;
                }
                let tm = tle.epoch + testvec[0] / 86400.0;

                // Test vectors assume WGS72 gravity model and AFSPC ops mode
                match sgp4_full(&mut tle, &[tm], GravConst::WGS72, OpsMode::AFSPC) {
                    Ok((pos, vel)) => {
                        for idx in 0..3 {
                            // Account for truncation in truth data
                            if testvec[idx + 4].abs() < 1.0e-4 {
                                maxvelerr = 1.0e-4;
                            }
                            if testvec[idx + 4].abs() < 1.0e-6 {
                                maxvelerr = 1.0e-2;
                            }
                            let poserr =
                                ((pos[idx] * 1.0e-3 - testvec[idx + 1]) / testvec[idx + 1]).abs();
                            let velerr =
                                ((vel[idx] * 1.0e-3 - testvec[idx + 4]) / testvec[idx + 4]).abs();
                            assert!(poserr < maxposerr);
                            assert!(velerr < maxvelerr);
                        }
                    }
                    // Note: some errors are part of test vectors
                    Err(e) => println!("SGP4 Error: \"{}\", is expected in testvecs", e.1),
                }
            }
        }
    }
}
