use super::astrotime::{AstroTime, Scale};
use std::f64::consts::PI;

use nalgebra as na;
type Vec3 = na::Vector3<f64>;
type Quat = na::UnitQuaternion<f64>;

use super::earth_orientation_params;
pub use crate::iau2000::qcirs2gcrs;

// Right-handed rotation of coordinate sytstem about x axis
// (left-handed rotation of vector)
#[inline]
pub(crate) fn qrotx(theta: f64) -> Quat {
    Quat::from_axis_angle(&Vec3::x_axis(), -theta)
}

// Right-handed rotation of coordinate sytstem about y axis
// (left-handed rotation of vector)
#[inline]
pub(crate) fn qroty(theta: f64) -> Quat {
    Quat::from_axis_angle(&Vec3::y_axis(), -theta)
}

// Right-handed rotation of coordinate sytstem about z axis
// (left-handed rotation of vector)
#[inline]
pub(crate) fn qrotz(theta: f64) -> Quat {
    Quat::from_axis_angle(&Vec3::z_axis(), -theta)
}

///
/// Greenwich Mean Sidereal Time
///
/// Vallado algorithm 15:
///
/// $\theta_{GMST} = 67310.5481 + (876600h + 8640184.812866) * t_{ut1} * (0.983104 + t_{ut1} * -6.2e-6)$
///
/// # Arguments
///
/// * tm: AstroTime object representing input time
///
/// # Returns
///
/// * gmst in radians
///
pub fn gmst(tm: &AstroTime) -> f64 {
    let tut1: f64 = (tm.to_mjd(Scale::UT1) - 51544.5) / 36525.0;
    let mut gmst: f64 = 67310.54841
        + tut1 * ((876600.0 * 3600.0 + 8640184.812866) + tut1 * (0.093104 + tut1 * -6.2e-6));

    gmst = (gmst % 86400.0) / 240.0 * PI / 180.0;
    return gmst;
}

///
/// Equation of Equinoxes
///
/// Vallado al
pub fn eqeq(tm: &AstroTime) -> f64 {
    let d: f64 = tm.to_mjd(Scale::TT) - 51544.5;
    let omega = PI / 180.0 * (125.04 - 0.052954 * d);
    let l = (280.47 + 0.98565 * d) * PI / 180.0;
    let epsilon = (23.4393 - 0.0000004 * d) * PI / 180.0;
    let d_psi = (-0.000319 * f64::sin(omega) - 0.000024 * f64::sin(2.0 * l)) * 15.0 * PI / 180.0;
    d_psi * f64::cos(epsilon)
}

pub fn gast(tm: &AstroTime) -> f64 {
    gmst(tm) + eqeq(tm)
}

///
/// Earth Rotation Angle
///
/// See
/// [IERS Technical Note 36, Chapter 5](https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36_043.pdf?__blob=publicationFile&v=1)
/// Equation 5.15
///
/// Arguments:
///
///   * tm: AstroTime struct representing input time
///
/// Returns:
///
///  * Earth rotation angle, in radians
///
pub fn earth_rotation_angle(tm: &AstroTime) -> f64 {
    let t = tm.to_jd(Scale::UT1);
    let f = t % 1.0;
    2.0 * PI * ((0.7790572732640 + f + 0.00273781191135448 * (t - 2451545.0)) % 1.0)
}

///
/// Rotation from International Terrestrial Reference Frame
/// to the Terrestrial Intermediate Reference System (TIRS)
///
/// # Arguments:
///  * tm: Time at which to compute rotation
///
/// # Return:
///
///  * Quaternion representing rotation from ITRF to TIRS
///
pub fn qitrf2tirs(tm: &AstroTime) -> Quat {
    const ASEC2RAD: f64 = PI / 180.0 / 3600.0;
    let eop = earth_orientation_params::get(tm).unwrap();
    let xp = eop[1] * ASEC2RAD;
    let yp = eop[2] * ASEC2RAD;
    let t_tt = (tm.to_mjd(Scale::TT) - 51544.5) / 36525.0;
    let sp = -47.0e-6 * ASEC2RAD * t_tt;
    qrotz(-sp) * qroty(xp) * qrotx(yp)
}

///
/// Rotation from True Equator Mean Equinox (TEME) frame
/// to International Terrestrial Reference Frame (TEME)
///
/// # Arguments
///
/// * tm: Time at which to compute rotation
///
/// # Returns
///
/// * Quaternion representing rotation from TEME to ITRF
///
/// # Notes
///
/// * The TEME frame is the default frame output by the
/// SGP4 propagator
/// * This is Equation 3-90 from Vallado
///
pub fn qteme2itrf(tm: &AstroTime) -> Quat {
    qitrf2tirs(tm).conjugate() * qrotz(gmst(tm))
}

///
/// Rotate from Mean Equinix of Date (MOD) coordinate frame
/// to Geocentric Celestrial Reference Frame
///
/// # Arguments
///
/// * tm: Time at which to compute rotation
///
/// # Returns
///
/// * Quaternion representing rotation from MOD to GCRF
///  
///
/// Equations 3-88 and 3-89 in Vallado
///
pub fn qmod2gcrs(tm: &AstroTime) -> Quat {
    const ASEC2RAD: f64 = PI / 180.0 / 3600.0;
    let tt = (tm.to_mjd(Scale::TT) - 51544.5) / 36525.0;

    let zeta = 2.650545
        + tt * (2306.083227
            + tt * (0.2988499 + tt * (0.01801828 + tt * (-0.000005971 + tt * -0.0000003173))));
    let z = -2.650545
        + tt * (2306.077181
            + tt * (1.0927348 + tt * (0.01826837 + tt * (-0.000028596 + tt * -0.0000002904))));
    let theta = tt
        * (2004.191903
            + tt * (-0.42949342 + tt * (-0.04182264 + tt * (-0.000007089 + tt * -0.0000001274))));

    return qrotz(zeta * ASEC2RAD) * qroty(-theta * ASEC2RAD) * qrotz(z * ASEC2RAD);
}

///
/// Approximate rotation from
/// Geocentric Celestrial Reference Frame to
/// International Terrestrial Reference Frame
///
/// This uses an approximation of the IAU-76/FK5 Reduction
/// See Vallado section 3.7.3
///
pub fn qgcrf2itrf_approx(tm: &AstroTime) -> Quat {
    // Neglecting polar motion
    let qitrf2tod_approx: Quat = qrotz(-gast(tm));

    qmod2gcrs(tm) * qtod2mod_approx(tm) * qitrf2tod_approx
}

/// Approximate rotation from
/// True of Date to Mean of Date
/// coordinate frame
///
/// See Vallado section 3.7.3
///
pub fn qtod2mod_approx(tm: &AstroTime) -> Quat {
    let d = tm.to_mjd(Scale::TT) - 51544.5;
    let t = d / 36525.0;

    const DEG2RAD: f64 = PI / 180.0;

    // Compute nutation rotation (accurate to ~ 1 arcsec)
    // This is where the approximation comes in
    let delta_psi = DEG2RAD
        * (-0.0048 * f64::sin((125.0 - 0.05295 * d) * DEG2RAD)
            - 0.0004 * f64::sin((200.9 + 1.97129 * d) * DEG2RAD));
    let delta_epsilon = DEG2RAD
        * (0.0026 * f64::cos((125.0 - 0.05295 * d) * DEG2RAD)
            + 0.0002 * f64::cos((200.9 + 1.97129 * d) * DEG2RAD));
    let epsilon_a = DEG2RAD
        * ((23.0 + 26.0 / 60.0 + 21.406 / 3600.0)
            + t * (-46.836769 / 3600.0
                + t * (-0.0001831 / 3600.0
                    + t * (0.00200340 / 3600.0
                        + t * (-5.76e-7 / 3600.0 + t * -4.34E-8 / 3600.0)))));
    let epsilon = epsilon_a + delta_epsilon;
    qrotx(-epsilon_a) * qrotz(delta_psi) * qrotx(epsilon)
}

///
/// Quaternion representing rotation from the
/// International Terrestrial Reference Frame (ITRF)
/// to the Geocentric Celestial Reference Frame (GCRF)
///
/// # Arguments
///
/// * tm: Time at which to compute rotation
///
/// # Returns
///
/// * Quaternion representing rotation from ITRF to GCRF
///
/// # Notes:
///
///  * Uses the full IAU2006 reduction, see
/// [IERS Technical Note 36, Chapter 5](https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36_043.pdf?__blob=publicationFile&v=1)
/// Equation 5.1
///
///  * This is **very** computationally expensive; for most
/// applications, the approximate rotation will work just fine
///
pub fn qitrf2gcrf(tm: &AstroTime) -> Quat {
    let w = qitrf2tirs(tm);
    let r = qtirs2cirs(tm);
    let q = qcirs2gcrs(tm);
    q * r * w
}

///
/// Quaternion representing rotation from the
/// Terrestrial Intermediate Reference System
/// to the Celestial Intermediate Reference System
///
/// A rotation about zhat by -Earth Rotation Angle
///
/// See [IERS Technical Note 36, Chapter 5](https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36_043.pdf?__blob=publicationFile&v=1)
/// Equation 5.5
///
#[inline]
pub fn qtirs2cirs(tm: &AstroTime) -> Quat {
    qrotz(-earth_rotation_angle(tm))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astrotime::{AstroTime, Scale};
    type Vec3 = na::Vector3<f64>;

    #[test]
    fn test_gmst() {
        // Vallado example 3-5
        let mut tm = AstroTime::from_datetime(1992, 8, 20, 12, 14, 0.0);
        // Spoof this as UT1 value
        let tdiff = tm.to_mjd(Scale::UT1) - tm.to_mjd(Scale::UTC);
        tm = tm - tdiff;
        // Convert to UT1
        let gmval = gmst(&tm) * 180.0 / PI;
        let truth = -207.4212121875;
        assert!(((gmval - truth) / truth).abs() < 1.0e-6)
    }

    #[test]
    fn test_gcrs2itrf() {
        // Example 3-14 from Vallado
        // With verification fo intermediate calculations
        // Input time
        let tm = &AstroTime::from_datetime(2004, 4, 6, 7, 51, 28.386009);
        // Input terrestrial location
        let pitrf = Vec3::new(-1033.4793830, 7901.2952754, 6380.3565958);
        let t_tt = (tm.to_jd(Scale::TT) - 2451545.0) / 36525.0;
        assert!((t_tt - 0.0426236319).abs() < 1.0e-8);

        let dut1 = (tm.to_mjd(Scale::UT1) - tm.to_mjd(Scale::UTC)) * 86400.0;
        // We linearly interpolate dut1, so this won't match exactly
        assert!((dut1 + 0.4399619).abs() < 0.01);
        let delta_at = (tm.to_mjd(Scale::TAI) - tm.to_mjd(Scale::UTC)) * 86400.0;
        assert!((delta_at - 32.0).abs() < 1.0e-7);
        let ptirs = qitrf2tirs(&tm) * pitrf;
        assert!((ptirs[0] + 1033.4750312).abs() < 1.0e-4);
        assert!((ptirs[1] - 7901.3055856).abs() < 1.0e-4);
        assert!((ptirs[2] - 6380.3445327).abs() < 1.0e-4);
        let era = earth_rotation_angle(&tm);
        assert!((era * 180.0 / PI - 312.7552829).abs() < 1.0e-5);
        let pcirs = qrotz(-era) * ptirs;
        assert!((pcirs[0] - 5100.0184047).abs() < 1e-3);
        assert!((pcirs[1] - 6122.7863648).abs() < 1e-3);
        assert!((pcirs[2] - 6380.3446237).abs() < 1e-3);
        let pgcrf = qcirs2gcrs(tm) * pcirs;
        assert!((pgcrf[0] - 5102.508959).abs() < 1e-3);
        assert!((pgcrf[1] - 6123.011403).abs() < 1e-3);
        assert!((pgcrf[2] - 6378.136925).abs() < 1e-3);
    }
}
