use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::{wrap_pyfunction, wrap_pymodule};

use std::path::PathBuf;

mod pyastrotime;
mod pycoordconversion;
mod pygravity;
mod pyitrfcoord;
mod pyjplephem;
mod pylpephem_moon;
mod pylpephem_sun;
mod pynrlmsise;
mod pyquaternion;
mod pysgp4;
mod pysolarsystem;
mod pytle;
mod pyuniv;

mod pyutils;

use pyastrotime::PyAstroTime;
use pycoordconversion as pycc;
use pyitrfcoord::PyITRFCoord;
use pyquaternion::Quaternion;
use pysolarsystem::SolarSystem;

/// JPL Ephemeris Sub-Module
#[pymodule]
fn jplephem(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pyjplephem::geocentric_pos, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pyjplephem::geocentric_state, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pyjplephem::barycentric_pos, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pyjplephem::barycentric_state, m)?)
        .unwrap();
    Ok(())
}

#[pymodule]
fn sun(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pylpephem_sun::pos_gcrf, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pylpephem_sun::pos_mod, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pylpephem_sun::rise_set, m)?)
        .unwrap();
    Ok(())
}

#[pymodule]
fn moon(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pylpephem_moon::pos_gcrf, m)?)
        .unwrap();
    Ok(())
}

/// Low-precision ephemerides
#[pymodule]
fn lpephem(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(sun))?;
    m.add_wrapped(wrap_pymodule!(moon))?;
    Ok(())
}

/// Return directory currently used to hold
/// necessary data files for the directory
///
/// e.g., Earth Orientation Parameters, gravity coefficients,
/// JPL Ephemeris, etc..
///
#[pyfunction]
fn datadir() -> PyResult<String> {
    Ok(String::from(
        crate::utils::datadir::get().unwrap().to_str().unwrap(),
    ))
}

///
/// Download data files needed for computation
///
/// Keyword Arguments:
///
///    overwrite:  <bool>  :: Download and overwrite files if they already exist
///          dir: <string> :: Target directory for files.  Uses existing
///                           data directory if not specifie
///
/// Files include:
///
///            EGM96.gfc :: EGM-96 Gravity Model Coefficients
///             JGM3.gfc :: JGM-3 Gravity Model Coefficients
///             JGM2.gfc :: JGM-2 Gravity Model Coefficients
///      ITU_GRACE16.gfc :: ITU Grace 16 Gravity Model Coefficients
///          tab5.2a.txt :: Coefficients for GCRS to GCRF conversion
///          tab5.2b.txt :: Coefficients for GCRS to GCRF conversion
///          tab5.2d.txt :: Coefficients for GCRS to GCRF conversion
///       sw19571001.txt :: Space weather data, updated daily
///     leap-seconds.txt :: Leap seconds (UTC vs TAI)
///      finals2000A.all :: Earth orientation parameters,  updated daily
/// linux_p1550p2650.440 :: JPL Ephemeris version 440 (~ 100 MB)
///
/// Note that files update daily will always be downloaded independed of
/// overwrite flag
///
#[pyfunction(kwds = "**")]
fn update_datafiles(kwds: Option<&PyDict>) -> PyResult<()> {
    let overwrite_files = match kwds {
        None => false,
        Some(u) => match u.get_item("overwrite") {
            Some(v) => v.extract::<bool>()?,
            None => false,
        },
    };
    let datadir = match kwds {
        None => None,
        Some(u) => match u.get_item("dir") {
            Some(v) => Some(PathBuf::from(v.extract::<String>()?)),
            None => None,
        },
    };

    match crate::utils::update_datafiles(datadir, overwrite_files) {
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e.to_string())),
        Ok(_) => Ok(()),
    }
}

#[pymodule]
fn util(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(datadir, m)?).unwrap();
    m.add_function(wrap_pyfunction!(update_datafiles, m)?)
        .unwrap();
    Ok(())
}

/// Coordinate Conversion Sub-Module
#[pymodule]
fn frametransform(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pycc::earth_rotation_angle, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pycc::gast, m)?).unwrap();
    m.add_function(wrap_pyfunction!(pycc::gmst, m)?).unwrap();
    m.add_function(wrap_pyfunction!(pycc::eqeq, m)?).unwrap();
    m.add_function(wrap_pyfunction!(pycc::qitrf2tirs, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pycc::qtirs2cirs, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pycc::qitrf2gcrf, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pycc::qteme2itrf, m)?)
        .unwrap();
    Ok(())
}

#[pymodule]
pub fn astro(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyAstroTime>()?;
    m.add_class::<pyastrotime::PyTimeScale>()?;
    m.add_class::<Quaternion>()?;
    m.add_function(wrap_pyfunction!(pysgp4::sgp4, m)?).unwrap();

    m.add_class::<pygravity::GravModel>()?;
    m.add_function(wrap_pyfunction!(pygravity::gravity, m)?)
        .unwrap();

    m.add_function(wrap_pyfunction!(pynrlmsise::nrlmsise00, m)?)
        .unwrap();

    m.add_class::<pyuniv::Univ>()?;
    m.add_class::<SolarSystem>()?;
    m.add_class::<pytle::PyTLE>()?;

    m.add_class::<PyITRFCoord>()?;

    m.add_wrapped(wrap_pymodule!(frametransform))?;
    m.add_wrapped(wrap_pymodule!(jplephem))?;
    m.add_wrapped(wrap_pymodule!(lpephem))?;

    m.add_wrapped(wrap_pymodule!(util))?;

    Ok(())
}
