use crate::univ;
use crate::AstroTime;
use crate::ITRFCoord;
use crate::TimeScale;

use crate::utils::{astroerr, AstroResult};

use nalgebra as na;

///
/// Sun position in the Geocentric Celestial Reference Frame (GCRF)
///
/// Algorithm 29 from Vallado for sun in Mean of Date (MOD), then rotated
/// from MOD to GCRF via Equations 3-88 and 3-89 in Vallado
///
/// Input:
///
///    time:  AstroTime object for which to compute position
///
/// Output:
///
///    nalagebra::Vector3<f64> representing sun position in GCRF frame
///    at given time.  Units are meters
///
#[inline]
pub fn pos_gcrf(time: &AstroTime) -> na::Vector3<f64> {
    crate::frametransform::qmod2gcrf(time) * pos_mod(time)
}

///
/// Sun position in the Mean-of-Date Frame
///
/// Algorithm 29 from Vallado for sun in Mean of Date (MOD)
///
/// Input:
///
///    time:  AstroTime object for which to compute position
///
/// Output:
///
///    nalagebra::Vector3<f64> representing sun position in MOD frame
///    at given time.  Units are meters
///
/// From Vallado: Valid with accuracy of .01 degrees from 1950 to 2050
///
pub fn pos_mod(time: &AstroTime) -> na::Vector3<f64> {
    let t: f64 = (time.to_jd(TimeScale::TDB) - 2451545.0) / 36525.0;
    #[allow(non_upper_case_globals)]
    const deg2rad: f64 = std::f64::consts::PI / 180.;

    // Mean longitude
    let lambda: f64 = 280.46 + 36000.77 * t;

    // mean anomaly
    #[allow(non_snake_case)]
    let M: f64 = deg2rad * (357.5277233 + 35999.05034 * t);

    // obliquity
    let epsilon: f64 = deg2rad * (23.439291 - 0.0130042 * t);

    // Ecliptic
    let lambda_ecliptic: f64 =
        deg2rad * (lambda + 1.914666471 * f64::sin(M) + 0.019994643 * f64::sin(2.0 * M));

    // Magnitude of sun vector
    let r: f64 =
        univ::AU * (1.000140612 - 0.016708617 * f64::cos(M) - 0.000139589 * f64::cos(2. * M));

    na::Vector3::<f64>::new(
        r * f64::cos(lambda_ecliptic),
        r * f64::sin(lambda_ecliptic) * f64::cos(epsilon),
        r * f64::sin(lambda_ecliptic) * f64::sin(epsilon),
    )
}

///
/// Sunrise and sunset times on the day given by input time
/// and at the given location.  
///
/// Time is at location, and should have hours, minutes, and seconds
/// set to zero
///
/// Vallado Algorithm 30
///
/// Inputs:
///
///      Time:   AstroTime representing date for which to compute
///              sunrise & sunset
///
///     coord:   ITRFCoord representing location for which to compute
///              sunrise & sunset
///
///     sigma:   Option<f64> Angle in degrees between noon & rise/set:
///              Common Values:
///                           "Standard": 90 deg, 50 arcmin (90.0+50.0/60.0)
///                     "Civil Twilight": 96 deg
///                  "Nautical Twilight": 102 deg
///              "Astronomical Twilight": 108 deg
///
///              If None is passed in, "Standard" is used (90.0 + 50.0/60.0)
///
/// Returns:
///
///    AstroResult<(sunrise: AstroTime, sunset: AstroTime)>
///
///
pub fn riseset(
    time: &AstroTime,
    coord: &ITRFCoord,
    osigma: Option<f64>,
) -> AstroResult<(AstroTime, AstroTime)> {
    use std::f64::consts::PI;
    let sigma = osigma.unwrap_or(90.0 + 50.0 / 60.0);
    let latitude: f64 = coord.latitude_deg();
    let longitude: f64 = coord.longitude_deg();

    let sind: fn(f64) -> f64 = |x: f64| (x * PI / 180.0).sin();
    let cosd: fn(f64) -> f64 = |x: f64| (x * PI / 180.0).cos();
    let tand: fn(f64) -> f64 = |x: f64| (x * PI / 180.0).tan();
    const RAD2DEG: f64 = 180.0 / PI;

    // Zero-hour GMST, equation 3-45 in Vallado
    let gmst0h = |t: f64| -> f64 {
        (100.4606184 + 36000.77005361 * t + 0.00038793 * t * t - 2.6E-8 * t * t * t) % 360.0
    };

    let jd0h = (time.to_jd(TimeScale::UTC) * 2.0).round() / 2.0;
    let criseset = |jdoffset: f64, lhafunc: fn(f64) -> f64| -> AstroResult<AstroTime> {
        let jd = time.to_jd(TimeScale::UTC) + jdoffset - longitude / 360.0;
        let t = (jd - 2451545.0) / 36525.0;

        let lambda_sun = 280.4606184 + 36000.77005361 * t;
        let msun = 357.5291092 + 35999.05034 * t;
        let lambda_ecliptic =
            lambda_sun + 1.914666471 * sind(msun) + 0.019994643 * sind(2.0 * msun);
        let epsilon = 23.439291 - 0.0130042 * t;

        let tanalpha_sun = cosd(epsilon) * tand(lambda_ecliptic);
        let sindelta_sun = sind(epsilon) * sind(lambda_ecliptic);
        let deltasun = f64::asin(sindelta_sun) * RAD2DEG;
        let alpha_sun = f64::atan(tanalpha_sun) * RAD2DEG;

        let coslha =
            (cosd(sigma) - sind(deltasun) * sind(latitude)) / (cosd(deltasun) * cosd(latitude));
        if coslha.abs() > 1.0 {
            return astroerr!(
                "Invalid position.  Sun doesn't rise/set on this day at \
                 this location (e.g., Alaska in summer)"
            );
        }
        let mut lha = f64::acos(coslha) * RAD2DEG;

        lha = lhafunc(lha);
        let gmst = gmst0h(t) % 360.0;
        let mut ret = (lha + alpha_sun - gmst) % 360.0;
        if ret < 0.0 {
            ret += 360.0;
        }
        Ok(AstroTime::from_jd(
            ret / 360.0 + jd0h - longitude / 360.0,
            TimeScale::UTC,
        ))
    };
    Ok((criseset(0.25, |x| 360.0 - x)?, criseset(0.75, |x| x)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sunpos_mod() {
        // Example 5-1 in Vallado
        let t0: AstroTime = AstroTime::from_date(2006, 4, 2);
        // Approximate this UTC as TDB to match example...
        let t = AstroTime::from_mjd(t0.to_mjd(TimeScale::UTC), TimeScale::TDB);

        let pos = pos_mod(&t);
        // Below value is from Vallado example
        let ref_pos = vec![146186212.0E3, 28788976.0E3, 12481064.0E3];
        for idx in 0..3 {
            let err = f64::abs(pos[idx] / ref_pos[idx] - 1.0);
            assert!(err < 1.0e-6);
        }
    }

    #[test]
    fn sunpos_gcrf() {
        // Example 5-1 in Vallado
        let t0: AstroTime = AstroTime::from_date(2006, 4, 2);
        // Approximate this UTC as TDB to match example...
        let t = AstroTime::from_mjd(t0.to_mjd(TimeScale::UTC), TimeScale::TDB);

        let pos = pos_gcrf(&t);
        // Below value is from Vallado example
        let ref_pos = vec![146259922.0E3, 28585947.0E3, 12397430.0E3];
        for idx in 0..3 {
            let err = f64::abs(pos[idx] / ref_pos[idx] - 1.0);
            // Less exact here because we are comparing to JPL ephemeris.
            // as described by Vallado
            assert!(err < 5e-4);
        }
    }

    #[test]
    fn sunriseset() {
        // Example 5-2 from Vallado
        let itrf = ITRFCoord::from_geodetic_deg(40.0, 0.0, 0.0);
        let tm = AstroTime::from_datetime(1996, 3, 23, 0, 0, 0.0);
        let (sunrise, sunset) = riseset(&tm, &itrf, None).unwrap();
        let (ryear, rmon, rday, rhour, rmin, rsec) = sunrise.to_datetime();
        assert!(ryear == 1996);
        assert!(rmon == 3);
        assert!(rday == 23);
        assert!(rhour == 5);
        assert!(rmin == 58);
        assert!((rsec / 21.97 - 1.0).abs() < 1.0e-3);
        let (syear, smon, sday, shour, smin, ssec) = sunset.to_datetime();
        assert!(syear == 1996);
        assert!(smon == 3);
        assert!(sday == 23);
        assert!(shour == 18);
        assert!(smin == 15);
        assert!((ssec / 17.76 - 1.0).abs() < 1.0e-3);

        // Check for error returned on 24-hour sunlight condition
        let itrf2 = ITRFCoord::from_geodetic_deg(85.0, 30.0, 0.0);
        let tm2 = AstroTime::from_date(2020, 6, 20);
        let r = riseset(&tm2, &itrf2, None);
        assert!(r.is_err());
    }
}
