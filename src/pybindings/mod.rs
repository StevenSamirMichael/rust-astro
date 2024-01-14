use pyo3::prelude::*;
use pyo3::{wrap_pyfunction, wrap_pymodule};

mod mod_utils;
mod pyastrotime;
mod pyconsts;
mod pydensity;
mod pyduration;
mod pyframetransform;
mod pygravity;
mod pyitrfcoord;
mod pyjplephem;
mod pylpephem_moon;
mod pylpephem_sun;
mod pynrlmsise;
mod pyquaternion;
mod pysatstate;
mod pysgp4;
mod pysolarsystem;
mod pytle;

mod pypropagate;
mod pypropsettings;
mod pysatproperties;

mod pyutils;

use pyastrotime::PyAstroTime;
use pyduration::PyDuration;
use pyframetransform as pyft;
use pyitrfcoord::PyITRFCoord;
use pyquaternion::Quaternion;
use pysolarsystem::SolarSystem;

use pypropsettings::PyPropSettings;
use pysatstate::PySatState;

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

/// Frame transform module: transform between varias coordinate frames
#[pymodule]
fn frametransform(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pyft::earth_rotation_angle, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pyft::gast, m)?).unwrap();
    m.add_function(wrap_pyfunction!(pyft::gmst, m)?).unwrap();
    m.add_function(wrap_pyfunction!(pyft::eqeq, m)?).unwrap();
    m.add_function(wrap_pyfunction!(pyft::qitrf2tirs, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pyft::qtirs2cirs, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pyft::qitrf2gcrf, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(pyft::qteme2itrf, m)?)
        .unwrap();
    Ok(())
}

#[pymodule]
fn satprop(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPropSettings>()?;
    m.add_class::<PySatState>()?;
    m.add_class::<pysatproperties::PySatProperties>()?;
    m.add_function(wrap_pyfunction!(pypropagate::propagate, m)?)
        .unwrap();

    Ok(())
}

#[pymodule] 
pub fn satkit(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyAstroTime>()?;
    m.add_class::<PyDuration>()?;
    m.add_class::<pyastrotime::PyTimeScale>()?;
    m.add_class::<Quaternion>()?;
    m.add_function(wrap_pyfunction!(pysgp4::sgp4, m)?).unwrap();

    m.add_class::<pygravity::GravModel>()?;
    m.add_function(wrap_pyfunction!(pygravity::gravity, m)?)
        .unwrap();

    m.add_function(wrap_pyfunction!(pynrlmsise::nrlmsise00, m)?)
        .unwrap();

    m.add_class::<pyconsts::Consts>()?;
    m.add_class::<SolarSystem>()?;
    m.add_class::<pytle::PyTLE>()?;

    m.add_class::<PyITRFCoord>()?;

    m.add_wrapped(wrap_pymodule!(frametransform))?;
    m.add_wrapped(wrap_pymodule!(jplephem))?;
    m.add_wrapped(wrap_pymodule!(sun))?;
    m.add_wrapped(wrap_pymodule!(moon))?;
    m.add_wrapped(wrap_pymodule!(satprop))?;

    m.add_wrapped(wrap_pymodule!(mod_utils::utils))?;
    m.add_wrapped(wrap_pymodule!(pydensity::density))?;

    Ok(())
}
