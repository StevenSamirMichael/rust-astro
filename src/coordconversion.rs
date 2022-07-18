use super::astrotime::{AstroTime, Scale};
use std::f64::consts::PI;

use nalgebra as na;
type Vec3 = na::Vector3<f64>;
type Quat = na::UnitQuaternion<f64>;

use super::earth_orientation_params;
use crate::iau2000::qcirs2gcrs;

#[inline]
pub fn qrotx(theta: f64) -> Quat {
    Quat::from_axis_angle(&Vec3::x_axis(), theta)
}

#[inline]
pub fn qroty(theta: f64) -> Quat {
    Quat::from_axis_angle(&Vec3::y_axis(), theta)
}

#[inline]
pub fn qrotz(theta: f64) -> Quat {
    Quat::from_axis_angle(&Vec3::z_axis(), theta)
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
/// See [IERS Technical Note 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html), Equation 5.15
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
    let t = tm.to_mjd(Scale::UT1);
    let f = t % 1.0;
    2.0 * PI * (0.7790572732640 + f + 0.00273781191135448 * t) % 1.0
}

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
/// Quaternion representing rotation from the
/// International Terrestrial Reference Frame (ITRF)
/// to the Geocentric Celestial Reference Frame (GCRS)
///
/// # Arguments
///
/// * tm: Time at which to compute transform
///
/// # Returns
///
/// * Quaternion representing rotation from ITRF to GCRS frame
///
/// Uses the IAU
pub fn qitrf2gcrs(tm: &AstroTime) -> Quat {
    let w = qitrf2tirs(tm);
    let r = qrotz(-earth_rotation_angle(tm));
    let q = qcirs2gcrs(tm);
    q * r * w
}

pub fn qgcrs2itrf(tm: &AstroTime) -> Quat {
    qitrf2gcrs(tm).conjugate()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astrotime::{AstroTime, Scale};

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
}
