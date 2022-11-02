use pyo3::prelude::*;
use pyo3::{wrap_pyfunction, wrap_pymodule};

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

#[pyfunction]
fn sum(a: i32, b: i32) -> PyResult<i32> {
    Ok(a + b)
}

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

#[pyfunction]
fn datadir() -> PyResult<String> {
    Ok(String::from(
        crate::utils::datadir::get().unwrap().to_str().unwrap(),
    ))
}

#[pymodule]
fn util(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(datadir, m)?).unwrap();
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
    m.add_wrapped(wrap_pyfunction!(sum))?;

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
