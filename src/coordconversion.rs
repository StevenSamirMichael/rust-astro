use super::astrotime::{AstroTime, Scale};
use std::f64::consts::PI;

use nalgebra as na;
type Vec3 = na::Vector3<f64>;
type Quat = na::UnitQuaternion<f64>;

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

pub fn gmst(tm: &AstroTime) -> f64 {
    let tut1: f64 = (tm.to_mjd(Scale::UT1) - 51544.5) / 36525.0;
    let mut gmst: f64 = 67310.54841
        + tut1 * ((876600.0 * 3600.0 + 8640184.812866) + tut1 * (0.093104 + tut1 * -6.2e-6));

    gmst = (gmst % 86400.0) / 240.0 * PI / 180.0;
    return gmst;
}

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
