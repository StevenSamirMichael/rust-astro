use super::pyquaternion::Quaternion;
use crate::astrotime::{AstroTime, TimeInput, TimeInputType};
use crate::frametransform::Quat;
use crate::utils::AstroResult;
use nalgebra as na;
use numpy as np;
use numpy::ndarray;
use numpy::{PyArray1, PyArray2};
use pyo3::prelude::*;

use crate::types::Vec3;

pub fn py_vec3_of_time_arr(
    cfunc: &dyn Fn(&AstroTime) -> Vec3,
    tmarr: &PyAny,
) -> PyResult<PyObject> {
    match tmarr.to_time_input() {
        TimeInput::Single(t) => {
            let v: Vec3 = cfunc(&t);
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                Ok(np::PyArray1::from_slice(py, &v.as_slice()).to_object(py))
            })
        }
        TimeInput::Array(tarr) => {
            let n = tarr.borrow().len();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                let out = np::PyArray2::<f64>::zeros(py, (n, 3), false);
                for idx in 0..n {
                    let v: Vec3 = cfunc(&tarr.borrow().get(idx).unwrap());
                    // I cannot figure out how to do this with a "safe" function,
                    // but... careful checking of dimensions above so this should
                    // never fail
                    unsafe {
                        std::ptr::copy_nonoverlapping(
                            v.as_ptr(),
                            out.as_raw_array_mut().as_mut_ptr().offset(idx as isize * 3),
                            3,
                        );
                    }
                }
                Ok(out.into_py(py))
            })
        }
        TimeInput::Error(_) => Err(pyo3::exceptions::PyTypeError::new_err(
            "Invalid input for time",
        )),
    }
}

pub fn py_vec3_of_time_result_arr(
    cfunc: &dyn Fn(&AstroTime) -> AstroResult<Vec3>,
    tmarr: &PyAny,
) -> PyResult<PyObject> {
    match tmarr.to_time_input() {
        TimeInput::Single(t) => match cfunc(&t) {
            Ok(v) => pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                Ok(np::PyArray1::from_slice(py, &v.as_slice()).to_object(py))
            }),
            Err(_) => Err(pyo3::exceptions::PyTypeError::new_err("Invalid time")),
        },
        TimeInput::Array(tarr) => {
            let n = tarr.borrow().len();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                let out = np::PyArray2::<f64>::zeros(py, (n, 3), false);
                for idx in 0..n {
                    match cfunc(&tarr.borrow().get(idx).unwrap()) {
                        Ok(v) => {
                            // I cannot figure out how to do this with a "safe" function,
                            // but... careful checking of dimensions above so this should
                            // never fail
                            unsafe {
                                std::ptr::copy_nonoverlapping(
                                    v.as_ptr(),
                                    out.as_raw_array_mut().as_mut_ptr().offset(idx as isize * 3),
                                    3,
                                );
                            }
                        }
                        Err(_) => {
                            return Err(pyo3::exceptions::PyTypeError::new_err("Invalid time"));
                        }
                    }
                }
                Ok(out.into_py(py))
            })
        }
        TimeInput::Error(_) => Err(pyo3::exceptions::PyTypeError::new_err(
            "Invalid input for time",
        )),
    }
}

pub fn py_func_of_time_arr<T: ToPyObject>(
    cfunc: fn(&AstroTime) -> T,
    tmarr: &PyAny,
) -> PyResult<PyObject> {
    match tmarr.to_time_input() {
        TimeInput::Single(v) => {
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> { Ok(cfunc(&v).to_object(py)) })
        }
        TimeInput::Array(arr) => {
            let tvec: Vec<T> = arr.borrow().iter().map(|x| cfunc(&x)).collect();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> { Ok(tvec.to_object(py)) })
        }
        TimeInput::Error(_) => Err(pyo3::exceptions::PyTypeError::new_err(
            "Invalid input for time",
        )),
    }
}

#[inline]
pub fn py_quat_from_time_arr(cfunc: fn(&AstroTime) -> Quat, tmarr: &PyAny) -> PyResult<PyObject> {
    match tmarr.to_time_input() {
        TimeInput::Single(v) => pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
            Ok(Quaternion { inner: cfunc(&v) }.into_py(py))
        }),
        TimeInput::Array(arr) => pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
            Ok(arr
                .borrow()
                .iter()
                .map(|x| -> Quaternion { Quaternion { inner: cfunc(&x) } })
                .collect::<Vec<Quaternion>>()
                .into_py(py))
        }),
        TimeInput::Error(_) => Err(pyo3::exceptions::PyTypeError::new_err(
            "Invalid input for time",
        )),
    }
}

#[inline]
pub fn tuple_func_of_time_arr<F>(cfunc: F, tmarr: &PyAny) -> PyResult<PyObject>
where
    F: Fn(&AstroTime) -> AstroResult<(na::Vector3<f64>, na::Vector3<f64>)>,
{
    match tmarr.to_time_input() {
        TimeInput::Single(v) => match cfunc(&v) {
            Ok(r) => pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                Ok((
                    PyArray1::from_slice(py, r.0.as_slice()),
                    PyArray1::from_slice(py, r.1.as_slice()),
                )
                    .to_object(py))
            }),
            Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e.to_string())),
        },
        TimeInput::Array(arr) => {
            let mut pout = ndarray::Array2::<f64>::zeros([arr.borrow().len(), 3]);
            let mut vout = ndarray::Array2::<f64>::zeros([arr.borrow().len(), 3]);
            for (i, tm) in arr.borrow().iter().enumerate() {
                match cfunc(tm) {
                    Ok(r) => {
                        pout.row_mut(i)
                            .assign(&ndarray::Array1::from_vec(vec![r.0[0], r.0[1], r.0[2]]));
                        vout.row_mut(i)
                            .assign(&ndarray::Array1::from_vec(vec![r.1[0], r.1[1], r.1[2]]));
                    }
                    Err(e) => return Err(pyo3::exceptions::PyRuntimeError::new_err(e.to_string())),
                }
            }
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                Ok((
                    PyArray2::from_array(py, &pout),
                    PyArray2::from_array(py, &vout),
                )
                    .to_object(py))
            })
        }
        TimeInput::Error(_) => Err(pyo3::exceptions::PyTypeError::new_err(
            "Invalid input for time",
        )),
    }
}
