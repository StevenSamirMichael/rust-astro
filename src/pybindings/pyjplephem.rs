use pyo3::prelude::*;

use super::pysolarsystem;
use super::pyutils::*;
use crate::astrotime::AstroTime;
use crate::jplephem::JPLEPHEM;
use crate::solarsystem::SolarSystem;
use crate::utils::AstroResult;
use nalgebra as na;

/// Return the position and velocity of the given body in
///  Geocentric coordinate system
///
/// # Inputs
///
///  * body - the solar system body for which to return position
///  * tm - The time at which to return position (can be list or array)
///
/// # Return
///
///   * Tuple with following elements:
///     * 3-vector of cartesian Geocentric position in meters
///       or Nx3 vector if input is list or array of times with length N
///
///     * 3-vector of cartesian Geocentric velocity in meters / second
///       or Nx3 vector if input is list or array of times with length N
///       Note: velocity is relative to Earth velocity
///
#[pyfunction]
pub fn geocentric_state(body: &pysolarsystem::SolarSystem, tm: &PyAny) -> PyResult<PyObject> {
    if JPLEPHEM.is_err() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(
            "Could not load ephemeris file",
        ));
    }

    let rbody: SolarSystem = body.into();
    let f = |tm: &AstroTime| -> AstroResult<(na::Vector3<f64>, na::Vector3<f64>)> {
        JPLEPHEM.as_ref().unwrap().geocentric_state(rbody, tm)
    };
    tuple_func_of_time_arr(f, tm)
}

/// Return the position & velocity the given body in the barycentric coordinate system
/// (origin is solar system barycenter)
///
/// # Inputs
///
///  * body - the solar system body for which to return position
///  * tm - The time at which to return position (can be list or array)
///
/// # Return
///
///  * Tuple with following values:
///    * 3-vector of cartesian Heliocentric position in meters
///      or Nx3 vector if input is list or array of times with length N
///
///    * 3-vector of cartesian Heliocentric velocity in meters / second
///      or Nx3 vector if input is list or array of times with length N
///
/// # Notes:
///  * Positions for all bodies are natively relative to solar system barycenter,
///    with exception of moon, which is computed in Geocentric system
///  * EMB (2) is the Earth-Moon barycenter
///  * The sun position is relative to the solar system barycenter
///    (it will be close to origin)
///
#[pyfunction]
pub fn barycentric_state(body: &pysolarsystem::SolarSystem, tm: &PyAny) -> PyResult<PyObject> {
    if JPLEPHEM.is_err() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(
            "Could not load ephemeris file",
        ));
    }

    let rbody: SolarSystem = body.into();
    let f = |tm: &AstroTime| -> AstroResult<(na::Vector3<f64>, na::Vector3<f64>)> {
        JPLEPHEM.as_ref().unwrap().barycentric_state(rbody, tm)
    };
    tuple_func_of_time_arr(f, tm)
}

/// Geocentric coordinate system
///
/// # Inputs
///
///  * body - the solar system body for which to return position
///  * tm - The time at which to return position (can be list or array)
///
/// # Return
///
///    3-vector of cartesian Geocentric position in meters
///    or Nx3 vector if input is list or array of times with length N
///
#[pyfunction]
pub fn geocentric_pos(body: &pysolarsystem::SolarSystem, tm: &PyAny) -> PyResult<PyObject> {
    if JPLEPHEM.is_err() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(
            "Could not load ephemeris file",
        ));
    }

    let rbody: SolarSystem = body.into();
    let f = |tm: &AstroTime| -> AstroResult<na::Vector3<f64>> {
        JPLEPHEM.as_ref().unwrap().geocentric_pos(rbody, tm)
    };
    py_vec3_of_time_result_arr(&f, tm)
}

/// Return the position of the given body in the Barycentric
/// coordinate system (origin is solarsystem barycenter)
///
/// # Inputs
///
///  * body - the solar system body for which to return position
///  * tm - The time at which to return position (can be list or array)
///
/// # Return
///
///    3-vector of cartesian Heliocentric position in meters
///    or Nx3 vector if input is list or array of times with length N
///
///
/// # Notes:
///  * Positions for all bodies are natively relative to solar system barycenter,
///    with exception of moon, which is computed in Geocentric system
///  * EMB (2) is the Earth-Moon barycenter
///  * The sun position is relative to the solar system barycenter
///    (it will be close to origin)
#[pyfunction]
pub fn barycentric_pos(body: &pysolarsystem::SolarSystem, tm: &PyAny) -> PyResult<PyObject> {
    if JPLEPHEM.is_err() {
        return Err(pyo3::exceptions::PyFileNotFoundError::new_err(
            "Could not load ephemeris file",
        ));
    }
    let rbody: SolarSystem = body.into();
    let f = |tm: &AstroTime| -> AstroResult<na::Vector3<f64>> {
        JPLEPHEM.as_ref().unwrap().barycentric_pos(rbody, tm)
    };
    py_vec3_of_time_result_arr(&f, tm)
}
