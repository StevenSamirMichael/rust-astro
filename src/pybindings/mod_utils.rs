use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::wrap_pyfunction;

use std::path::PathBuf;

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
#[pyfunction]
#[pyo3(signature=(**kwds))]
fn update_datafiles(kwds: Option<&PyDict>) -> PyResult<()> {
    let overwrite_files = match kwds {
        None => false,
        Some(u) => match u.get_item("overwrite")? {
            Some(v) => v.extract::<bool>()?,
            None => false,
        },
    };
    let datadir = match kwds {
        None => None,
        Some(u) => match u.get_item("dir")? {
            Some(v) => Some(PathBuf::from(v.extract::<String>()?)),
            None => None,
        },
    };

    match crate::utils::update_datafiles(datadir, overwrite_files) {
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e.to_string())),
        Ok(_) => Ok(()),
    }
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

/// Git hash of compiled library
#[pyfunction]
fn githash() -> PyResult<String> {
    Ok(String::from(crate::utils::githash()))
}

/// Build date of compiled library
#[pyfunction]
fn build_date() -> PyResult<String> {
    Ok(String::from(crate::utils::build_date()))
}

/// Astro utility functions
#[pymodule]
pub fn utils(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(datadir, m)?).unwrap();
    m.add_function(wrap_pyfunction!(update_datafiles, m)?)
        .unwrap();
    m.add_function(wrap_pyfunction!(githash, m)?).unwrap();
    m.add_function(wrap_pyfunction!(build_date, m)?).unwrap();
    Ok(())
}
