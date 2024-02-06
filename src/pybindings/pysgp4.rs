use pyo3::prelude::*;

use super::pyastrotime::ToTimeVec;
use super::pytle::PyTLE;
use crate::sgp4 as psgp4;
use numpy::PyArray1;

#[pyfunction]
pub fn sgp4(tle_input: &PyAny, tm: &PyAny) -> PyResult<PyObject> {
    if tle_input.is_instance_of::<PyTLE>() {
        let mut tle: PyRefMut<PyTLE> = tle_input.extract()?;
        match psgp4::sgp4(&mut tle.inner, tm.to_time_vec()?.as_slice()) {
            Ok((r, v)) => pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
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
                )
                    .to_object(py))
            }),
            Err(e) => {
                let estr = format!("Error running sgp4: {}", e.1);
                Err(pyo3::exceptions::PyRuntimeError::new_err(estr))
            }
        }
    } else {
        Err(pyo3::exceptions::PyRuntimeError::new_err(
            "Invalid input type for argument 1",
        ))
    }
}
