use super::pyastrotime::PyAstroTime;
use super::pypropsettings::PyPropSettings;

use nalgebra as na;
use numpy as np;

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::orbitprop::{PropSettings, SatState};

#[pyclass(name = "SatState")]
#[derive(Clone, Debug)]
pub struct PySatState {
    inner: SatState,
}

#[pymethods]
impl PySatState {
    #[new]
    fn py_new(
        tm: &PyAstroTime,
        pos: &np::PyArray1<f64>,
        vel: &np::PyArray1<f64>,
    ) -> PyResult<Self> {
        if pos.len() != 3 || vel.len() != 3 {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Position and velocity must be 1-d numpy arrays with length 3",
            ));
        }

        Ok(PySatState {
            inner: SatState::from_pv(
                &tm.inner,
                &na::Vector3::<f64>::from_row_slice(unsafe { pos.as_slice().unwrap() }),
                &na::Vector3::<f64>::from_row_slice(unsafe { vel.as_slice().unwrap() }),
            ),
        })
    }

    #[pyo3(signature=(tm, **kwargs))]
    fn propagate(&self, tm: &PyAstroTime, kwargs: Option<&PyDict>) -> PyResult<Self> {
        let propsettings: Option<PropSettings> = match kwargs.is_some() {
            true => {
                let kw = kwargs.unwrap();
                match kw.get_item("propsettings")? {
                    None => None,
                    Some(v) => Some(v.extract::<PyPropSettings>()?.inner),
                }
            }
            false => None,
        };

        match self.inner.propagate(&tm.inner, propsettings.as_ref()) {
            Ok(s) => Ok(PySatState { inner: s }),
            Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(format!(
                "Error propagating state: {}",
                e.to_string()
            ))),
        }
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner.to_string())
    }
}
