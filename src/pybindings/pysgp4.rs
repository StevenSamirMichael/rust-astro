use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use super::pyastrotime::ToTimeVec;
use super::pytle::PyTLE;
use crate::sgp4 as psgp4;
use numpy::PyArray1;

#[pyfunction]
#[pyo3(signature=(tle, time, **kwds))]
pub fn sgp4(tle: &PyAny, time: &PyAny, kwds: Option<&PyDict>) -> PyResult<PyObject> {
    let mut output_err = false;
    if kwds.is_some() {
        let kw = kwds.unwrap();
        match kw.get_item("errflag").unwrap() {
            Some(v) => output_err = v.extract::<bool>()?,
            None => {}
        }
    }
    if tle.is_instance_of::<PyTLE>() {
        let mut stle: PyRefMut<PyTLE> = tle.extract()?;
        match psgp4::sgp4(&mut stle.inner, time.to_time_vec()?.as_slice()) {
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
                if output_err == true {
                    pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                        Ok((e.0, e.1).to_object(py))
                    })
                } else {
                    let estr = format!("Error running sgp4: {}", e.1);
                    Err(pyo3::exceptions::PyRuntimeError::new_err(estr))
                }
            }
        }
    } else if tle.is_instance_of::<PyList>() {
        let mut tles = tle.extract::<Vec<PyRefMut<PyTLE>>>()?;
        let tmarray = time.to_time_vec()?;
        let results: Vec<psgp4::SGP4Result> = tles
            .iter_mut()
            .map(|tle| psgp4::sgp4(&mut tle.inner, tmarray.as_slice()))
            .collect();
        pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
            let n = tles.len() * tmarray.len() * 3;
            let mut parr: &PyArray1<f64> = PyArray1::zeros(py, [n], false);
            let mut varr: &PyArray1<f64> = PyArray1::zeros(py, [n], false);
            results.iter().enumerate().for_each(|(idx, r)| {
                let x = 3;
            });
            Ok((parr, varr).to_object(py))
        })
    } else {
        Err(pyo3::exceptions::PyRuntimeError::new_err(
            "Invalid input type for argument 1",
        ))
    }
}
