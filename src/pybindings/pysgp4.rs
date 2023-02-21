use pyo3::prelude::*;

use super::pyastrotime::ToTimeVec;
use super::pytle::PyTLE;
use crate::sgp4 as psgp4;
use numpy::PyArray1;

#[pyfunction]
pub fn sgp4(tle: &mut PyTLE, tm: &PyAny) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
    match psgp4::sgp4(tle.into(), tm.to_time_vec()?.as_slice()) {
        Ok((r, v)) => pyo3::Python::with_gil(|py| -> PyResult<(Py<PyAny>, Py<PyAny>)> {
            let mut dims = vec![r.len()];
            if r.nrows() > 1 && r.ncols() > 1 {
                dims = vec![r.ncols(), r.nrows()];
            }

            // Note: this is a little confusing: ndarray uses
            // row major, nalgebra and numpy use column major,
            // hence the switch
            Ok((
                PyArray1::from_slice(py, r.data.as_slice())
                    .reshape(dims.clone())
                    .unwrap()
                    .to_object(py),
                PyArray1::from_slice(py, v.data.as_slice())
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
