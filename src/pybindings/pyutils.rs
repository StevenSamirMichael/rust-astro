use super::pyastrotime::ToTimeVec;
use super::pyquaternion::Quaternion;

use crate::astrotime::AstroTime;
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
    let tm = tmarr.to_time_vec()?;
    match tm.len() {
        1 => {
            let v: Vec3 = cfunc(&tm[0]);
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                Ok(np::PyArray1::from_slice(py, &v.as_slice()).to_object(py))
            })
        }
        _ => {
            let n = tm.len();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                let out = np::PyArray2::<f64>::zeros(py, (n, 3), false);
                for idx in 0..n {
                    let v: Vec3 = cfunc(&tm[idx]);
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
    }
}

pub fn py_vec3_of_time_result_arr(
    cfunc: &dyn Fn(&AstroTime) -> AstroResult<Vec3>,
    tmarr: &PyAny,
) -> PyResult<PyObject> {
    let tm = tmarr.to_time_vec()?;

    match tm.len() {
        1 => match cfunc(&tm[0]) {
            Ok(v) => pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                Ok(np::PyArray1::from_slice(py, &v.as_slice()).to_object(py))
            }),
            Err(_) => Err(pyo3::exceptions::PyTypeError::new_err("Invalid time")),
        },
        _ => {
            let n = tm.len();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                let out = np::PyArray2::<f64>::zeros(py, (n, 3), false);
                for idx in 0..n {
                    match cfunc(&tm[idx]) {
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
    }
}

pub fn py_func_of_time_arr<T: ToPyObject>(
    cfunc: fn(&AstroTime) -> T,
    tmarr: &PyAny,
) -> PyResult<PyObject> {
    let tm = tmarr.to_time_vec()?;

    match tm.len() {
        1 => pyo3::Python::with_gil(|py| -> PyResult<PyObject> { Ok(cfunc(&tm[0]).to_object(py)) }),
        _ => {
            let tvec: Vec<T> = tm.iter().map(|x| cfunc(&x)).collect();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> { Ok(tvec.to_object(py)) })
        }
    }
}

#[inline]
pub fn py_quat_from_time_arr(cfunc: fn(&AstroTime) -> Quat, tmarr: &PyAny) -> PyResult<PyObject> {
    let tm = tmarr.to_time_vec()?;
    match tm.len() {
        1 => pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
            Ok(Quaternion {
                inner: cfunc(&tm[0]),
            }
            .into_py(py))
        }),
        _ => pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
            Ok(tm
                .iter()
                .map(|x| -> Quaternion { Quaternion { inner: cfunc(&x) } })
                .collect::<Vec<Quaternion>>()
                .into_py(py))
        }),
    }
}

#[inline]
pub fn tuple_func_of_time_arr<F>(cfunc: F, tmarr: &PyAny) -> PyResult<PyObject>
where
    F: Fn(&AstroTime) -> AstroResult<(na::Vector3<f64>, na::Vector3<f64>)>,
{
    let tm = tmarr.to_time_vec()?;
    match tm.len() {
        1 => match cfunc(&tm[0]) {
            Ok(r) => pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                Ok((
                    PyArray1::from_slice(py, r.0.as_slice()),
                    PyArray1::from_slice(py, r.1.as_slice()),
                )
                    .to_object(py))
            }),
            Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e.to_string())),
        },
        _ => {
            let mut pout = ndarray::Array2::<f64>::zeros([tm.len(), 3]);
            let mut vout = ndarray::Array2::<f64>::zeros([tm.len(), 3]);
            for (i, tm) in tm.iter().enumerate() {
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
    }
}
