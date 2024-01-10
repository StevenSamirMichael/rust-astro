use pyo3::prelude::*;

use crate::consts as cconsts;

#[pyclass(name = "consts")]
pub struct Consts {}

#[pymethods]
#[allow(non_upper_case_globals)]
impl Consts {
    #[classattr]
    ///  WGS-84 semiparameter, in meters
    const wgs84_a: f64 = cconsts::WGS84_A;
    #[classattr]
    ///  WGS-84 flattening
    const wgs84_f: f64 = cconsts::WGS84_F;
    #[classattr]
    /// WGS-84 Earth radius, meters
    const earth_radius: f64 = cconsts::EARTH_RADIUS;

    ///  Gravitational parameter of Earth, m^3/s^2 */
    #[classattr]
    const mu_earth: f64 = cconsts::MU_EARTH;
    #[classattr]
    const mu_moon: f64 = cconsts::MU_MOON;
    #[classattr]
    const mu_sun: f64 = cconsts::MU_SUN;
    #[classattr]
    const GM: f64 = cconsts::GM;
    #[classattr]
    const omega_earth: f64 = cconsts::OMEGA_EARTH;
    #[classattr]
    const c: f64 = cconsts::C;
    #[classattr]
    const au: f64 = cconsts::AU;
    #[classattr]
    const sun_radius: f64 = cconsts::SUN_RADIUS;
    #[classattr]
    const moon_radius: f64 = cconsts::MOON_RADIUS;
    #[classattr]
    const earth_moon_mass_ratio: f64 = cconsts::EARTH_MOON_MASS_RATIO;
    #[classattr]
    const geo_r: f64 = cconsts::GEO_R;
    #[classattr]
    const jgm3_mu: f64 = cconsts::JGM3_MU;
    #[classattr]
    const jgm3_a: f64 = cconsts::JGM3_A;
    #[classattr]
    const jgm3_j2: f64 = cconsts::JGM3_J2;
}
