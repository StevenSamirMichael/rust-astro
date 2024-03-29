use pyo3::prelude::*;

use crate::gravity::{accel, GravityModel};

use super::pyitrfcoord::PyITRFCoord;
use crate::itrfcoord::ITRFCoord;
use nalgebra as na;
use numpy as np;

use pyo3::types::PyDict;

///
/// Gravity model enumeration
///
/// For details of models, see:
/// http://icgem.gfz-potsdam.de/tom_longtime
///
#[allow(non_camel_case_types)]
#[pyclass(name = "gravmodel")]
#[derive(Clone)]
pub enum GravModel {
    jgm3 = GravityModel::JGM3 as isize,
    jgm2 = GravityModel::JGM2 as isize,
    egm96 = GravityModel::EGM96 as isize,
    itugrace16 = GravityModel::ITUGrace16 as isize,
}

impl From<GravModel> for GravityModel {
    fn from(g: GravModel) -> GravityModel {
        match g {
            GravModel::jgm3 => GravityModel::JGM3,
            GravModel::jgm2 => GravityModel::JGM2,
            GravModel::egm96 => GravityModel::EGM96,
            GravModel::itugrace16 => GravityModel::ITUGrace16,
        }
    }
}

impl IntoPy<PyObject> for &GravityModel {
    fn into_py(self, py: Python<'_>) -> PyObject {
        let g: GravModel = match self {
            GravityModel::JGM3 => GravModel::jgm3,
            GravityModel::JGM2 => GravModel::jgm2,
            GravityModel::EGM96 => GravModel::egm96,
            GravityModel::ITUGrace16 => GravModel::itugrace16,
        };
        g.into_py(py)
    }
}

///
/// gravity(pos, /)
/// --
///
/// Return acceleration due to Earth gravity at the input position. The
/// acceleration does not include the centrifugal force, and is output
/// in m/s^2 in the International Terrestrial Reference Frame (ITRF)
///
/// Inputs:
///
///       pos:   Position as ITRF coordinate (astro.itrfcoord) or numpy
///              3-vector representing ITRF position in meters
///
/// Kwargs:
///     
///     model:   The gravity model to use.  Options are:
///                   astro.GravModel.jgm3
///                   astro.GravModel.jgm2
///                   astro.GravModel.egm96
///                   astro.GravModel.itugrace16
///
///               Default is astro.GravModel.jgm3
///
///               For details of models, see:
///               http://icgem.gfz-potsdam.de/tom_longtime
///
///     order:    The order of the gravity model to use.
///               Default is 6, maximum is 16
///
///
///               For details of calculation, see Chapter 3.2 of:
///               "Satellite Orbits: Models, Methods, Applications",
///               O. Montenbruck and B. Gill, Springer, 2012.
///
#[pyfunction]
#[pyo3(signature=(pos, **kwds))]
pub fn gravity(pos: &PyAny, kwds: Option<&PyDict>) -> PyResult<PyObject> {
    let mut order: usize = 6;
    let mut model: GravModel = GravModel::jgm3;
    if kwds.is_some() {
        let kw = kwds.unwrap();
        if kw.contains("model").unwrap() {
            match kw.get_item("model")?.unwrap().extract::<GravModel>() {
                Ok(v) => model = v,
                Err(e) => return Err(e),
            }
        }
        if kw.contains("order").unwrap() {
            match kw.get_item("order")?.unwrap().extract::<usize>() {
                Ok(v) => order = v,
                Err(e) => return Err(e),
            }
        }
    }

    if pos.is_instance_of::<PyITRFCoord>() {
        let pyitrf: PyRef<PyITRFCoord> = pos.extract()?;
        let itrf: ITRFCoord = pyitrf.inner.into();
        let v = accel(&itrf.itrf, order, model.into());
        pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
            let vpy: &np::PyArray1<f64> = np::PyArray1::<f64>::from_slice(py, v.as_slice());
            Ok(vpy.into_py(py))
        })
    } else if pos.is_instance_of::<np::PyArray1<f64>>() {
        let vpy = pos.extract::<np::PyReadonlyArray1<f64>>().unwrap();
        if vpy.len() != 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "input must have 3 elements",
            ));
        }
        let v: na::Vector3<f64> = na::Vector3::<f64>::from_row_slice(vpy.as_slice().unwrap());
        let a = accel(&v, order, model.into());
        pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
            let vpy = np::PyArray1::<f64>::from_slice(py, a.as_slice());
            Ok(vpy.into_py(py))
        })
    } else {
        Err(pyo3::exceptions::PyTypeError::new_err(
            "Input must be 3-elemnent numpy or itrfcoord",
        ))
    }
}
