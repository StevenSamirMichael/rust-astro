use pyo3::prelude::*;

use crate::univ as cuniv;

#[pyclass(name = "consts")]
pub struct Consts {}

#[pymethods]
#[allow(non_upper_case_globals)]
impl Consts {
    #[classattr]
    ///  WGS-84 semiparameter, in meters
    const wgs84_a: f64 = cuniv::WGS84_A;
    #[classattr]
    ///  WGS-84 flattening
    const wgs84_f: f64 = cuniv::WGS84_F;
    #[classattr]
    /// WGS-84 Earth radius, meters
    const earth_radius: f64 = cuniv::EARTH_RADIUS;

    ///  Gravitational parameter of Earth, m^3/s^2 */
    #[classattr]
    const mu_earth: f64 = cuniv::MU_EARTH;
    #[classattr]
    const mu_moon: f64 = cuniv::MU_MOON;
    #[classattr]
    const mu_sun: f64 = cuniv::MU_SUN;
    #[classattr]
    const GM: f64 = cuniv::GM;
    #[classattr]
    const omega_earth: f64 = cuniv::OMEGA_EARTH;
    #[classattr]
    const c: f64 = cuniv::C;
    #[classattr]
    const au: f64 = cuniv::AU;
    #[classattr]
    const sun_radius: f64 = cuniv::SUN_RADIUS;
    #[classattr]
    const moon_radius: f64 = cuniv::MOON_RADIUS;
    #[classattr]
    const earth_moon_mass_ratio: f64 = cuniv::EARTH_MOON_MASS_RATIO;
    #[classattr]
    const geo_r: f64 = cuniv::GEO_R;
    #[classattr]
    const jgm3_mu: f64 = cuniv::JGM3_MU;
    #[classattr]
    const jgm3_a: f64 = cuniv::JGM3_A;
    #[classattr]
    const jgm3_j2: f64 = cuniv::JGM3_J2;
}
