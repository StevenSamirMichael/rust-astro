use super::sgp4_lowlevel::sgp4_lowlevel;
use super::sgp4init::sgp4init;

use crate::astrotime::AstroTime;
use crate::astrotime::Scale;
use crate::tle::TLE;

type Vec3 = nalgebra::Vector3<f64>;

type SGP4State = (Vec3, Vec3);
type SGP4Result = Result<SGP4State, i32>;

use std::f64::consts::PI;

pub fn sgp4(tle: &mut TLE, tm: AstroTime) -> SGP4Result {
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

        tle.satrec = Some(sgp4init(
            &"wgs84",
            'i',
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
        ));
    }
    let s = tle.satrec.as_mut().unwrap();

    let mut r: [f64; 3] = [0.0, 0.0, 0.0];
    let mut v: [f64; 3] = [0.0, 0.0, 0.0];
    let tsince = (tm - tle.epoch) * 1440.0;
    sgp4_lowlevel(s, tsince, &mut r, &mut v);

    if s.error != 0 {
        return Err(s.error as i32);
    }
    Ok((
        Vec3::new(r[0], r[1], r[2]) * 1.0e3,
        Vec3::new(v[0], v[1], v[2]) * 1.0e3,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tle::TLE;

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
        println!("tle = {:?}", tle);
        let tm = tle.epoch;
        println!("tm = {}", tm);

        match sgp4(&mut tle, tm) {
            Ok((pos, vel)) => {
                println!("pos = {}", pos);
                println!("vel = {}", vel);
            }
            Err(e) => {
                panic!("Error running sgp4: {}", e);
            }
        }
    }
}
