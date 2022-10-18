use pyo3::prelude::*;

use super::pytle::PyTLE;
use crate::sgp4 as psgp4;
use numpy::PyArray1;

#[pyfunction]
pub fn sgp4(tle: &PyTLE, tm: &PyAny) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
    match psgp4::sgp4(tle.into(), &tm) {
        Ok(v) => pyo3::Python::with_gil(|py| -> PyResult<(Py<PyAny>, Py<PyAny>)> {
            let mut dims = vec![v.0.len()];
            if v.0.nrows() > 1 && v.0.ncols() > 1 {
                dims = vec![v.0.ncols(), v.0.nrows()];
            }

            // Note: this is a little confusing: ndarray uses
            // row major, nalgebra and numpy use column major,
            // hence the switch
            Ok((
                PyArray1::from_slice(py, v.0.data.as_slice())
                    .reshape(dims.clone())
                    .unwrap()
                    .to_object(py),
                PyArray1::from_slice(py, v.1.data.as_slice())
                    .reshape(dims)
                    .unwrap()
                    .to_object(py),
            ))
        }),
        Err(e) => {
            let estr = format!("Error running sgp4: {}", e.1);
            Err(pyo3::exceptions::PyRuntimeError::new_err(estr))
        }
    }
}
