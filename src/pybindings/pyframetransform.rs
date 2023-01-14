use super::pyutils::*;
use crate::frametransform as ft;
use pyo3::prelude::*;

///
/// Greenwich Mean Sidereal Time
///
/// Vallado algorithm 15:
///
/// GMST = 67310.5481 + (876600h + 8640184.812866) * tᵤₜ₁ * (0.983104 + tᵤₜ₁ * −6.2e−6)
///
///
/// # Arguments
///
/// * tm: AstroTime object representing input time
///
/// # Returns
///
/// * Greenwich Mean Sideral Time,  radians
///
#[pyfunction]
pub fn gmst(tm: &PyAny) -> PyResult<PyObject> {
    py_func_of_time_arr(ft::gmst, tm)
}

///
/// Equation of Equinoxes
///
#[pyfunction]
pub fn eqeq(tm: &PyAny) -> PyResult<PyObject> {
    py_func_of_time_arr(ft::eqeq, tm)
}

///
/// Greenwich apparant sidereal time, radians
///
/// # Arguments:
///
///   * tm: astro.time struct representing input time
///
/// # Returns:
///
///  * Greenwich apparant sidereal time, radians
///
#[pyfunction]
pub fn gast(tm: &PyAny) -> PyResult<PyObject> {
    py_func_of_time_arr(ft::gast, tm)
}

///
/// Earth Rotation Angle
///
/// See
/// [IERS Technical Note 36, Chapter 5](https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36_043.pdf?__blob=publicationFile&v=1)
/// Equation 5.15
///
/// # Arguments:
///
///   * tm: AstroTime struct representing input time
///
/// # Returns:
///
///  * Earth rotation angle, in radians
///
/// # Calculation Details
///
/// * Let t be UT1 Julian date
/// * let f be fractional component of t (fraction of day)
/// * ERA = 2𝜋 ((0.7790572732640 + f + 0.00273781191135448 * (t − 2451545.0))
///
///
#[pyfunction]
pub fn earth_rotation_angle(tm: &PyAny) -> PyResult<PyObject> {
    py_func_of_time_arr(ft::earth_rotation_angle, tm)
}

///
/// Rotation from International Terrestrial Reference Frame
/// (ITRF) to the Terrestrial Intermediate Reference System (TIRS)
///
/// # Arguments:
///
///   * tm: astro.time struct representing input time
///
/// # Returns:
///
///  * Quaternion representing rotation from ITRF to TIRS
///
#[pyfunction]
pub fn qitrf2tirs(tm: &PyAny) -> PyResult<PyObject> {
    py_quat_from_time_arr(ft::qitrf2tirs, tm)
}

#[pyfunction]
pub fn qtirs2cirs(tm: &PyAny) -> PyResult<PyObject> {
    py_quat_from_time_arr(ft::qtirs2cirs, tm)
}

///
/// Quaternion representing rotation from the
/// International Terrestrial Reference Frame (ITRF)
/// to the Geocentric Celestial Reference Frame (GCRF)
///
/// Uses full IAU2006 Reduction
/// See IERS Technical Note 36, Chapter 5
///
/// Note: Very computationally expensive
///
/// # Arguments:
///
///   * tm: astro.time struct representing input time
///
/// # Returns:
///
///  * Quaternion representing rotation from ITRF to GCRF
///
#[pyfunction]
pub fn qitrf2gcrf(tm: &PyAny) -> PyResult<PyObject> {
    py_quat_from_time_arr(ft::qitrf2gcrf, tm)
}

///
/// Rotation from True Equator Mean Equinox (TEME) frame
/// to International Terrestrial Reference Frame (ITRF)
///
/// Note: TEME is output frame of SGP4 propagator
///
/// This is Equation 3-90 in Vallado
///
/// # Arguments:
///
///   * tm: astro.time struct representing input time
///
/// # Returns:
///
///  * Quaternion representing rotation from TEME to ITRF
///
#[pyfunction]
pub fn qteme2itrf(tm: &PyAny) -> PyResult<PyObject> {
    py_quat_from_time_arr(ft::qteme2itrf, tm)
}
