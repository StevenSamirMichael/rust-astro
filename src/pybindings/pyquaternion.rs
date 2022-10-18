use nalgebra as na;
use numpy as np;
use numpy::ndarray::array;
use numpy::ToPyArray;
use pyo3::prelude::*;

type Quat = na::UnitQuaternion<f64>;
type Vec3 = na::Vector3<f64>;

#[pyclass]
#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Quaternion {
    pub inner: Quat,
}

#[pymethods]
impl Quaternion {
    #[new]
    fn py_new() -> PyResult<Self> {
        Ok(Quaternion {
            inner: Quat::from_axis_angle(&Vec3::x_axis(), 0.0),
        })
    }

    #[staticmethod]
    fn rotx(theta_rad: f64) -> PyResult<Self> {
        Ok(Quaternion {
            inner: Quat::from_axis_angle(&Vec3::x_axis(), theta_rad),
        })
    }

    #[staticmethod]
    fn roty(theta_rad: f64) -> PyResult<Self> {
        Ok(Quaternion {
            inner: Quat::from_axis_angle(&Vec3::y_axis(), theta_rad),
        })
    }

    #[staticmethod]
    fn rotz(theta_rad: f64) -> PyResult<Self> {
        Ok(Quaternion {
            inner: Quat::from_axis_angle(&Vec3::z_axis(), theta_rad),
        })
    }

    #[staticmethod]
    fn from_axis_angle(axis: np::PyReadonlyArray1<f64>, angle: f64) -> PyResult<Self> {
        let v: Vec3 = Vec3::new(axis.as_array()[0], axis.as_array()[1], axis.as_array()[2]);
        let u = na::UnitVector3::try_new(v, 1.0e-9);
        match u {
            Some(ax) => Ok(Quaternion {
                inner: Quat::from_axis_angle(&ax, angle),
            }),
            None => {
                let err = pyo3::exceptions::PyArithmeticError::new_err("Axis norm is 0");
                Err(err)
            }
        }
    }

    fn __str__(&self) -> PyResult<String> {
        let ax: na::Unit<Vec3> = match self.inner.axis() {
            Some(v) => v,
            None => na::Unit::new_normalize(Vec3::new(1.0, 0.0, 0.0)),
        };
        let angle = self.inner.angle();
        Ok(format!(
            "Quaternion(Axis = [{:6.4}, {:6.4}, {:6.4}], Angle = {:6.4} rad)",
            ax[0], ax[1], ax[2], angle
        ))
    }

    fn __repr__(&self) -> PyResult<String> {
        self.__str__()
    }

    #[getter]
    fn angle(&self) -> PyResult<f64> {
        Ok(self.inner.angle())
    }

    #[getter]
    fn axis(&self) -> PyResult<PyObject> {
        let a = match self.inner.axis() {
            Some(ax) => ax,
            None => Vec3::x_axis(),
        };
        pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
            Ok(array![a[0], a[1], a[2]].to_pyarray(py).to_object(py))
        })
    }

    fn conj(&self) -> PyResult<Quaternion> {
        Ok(Quaternion {
            inner: self.inner.conjugate(),
        })
    }

    fn conjugate(&self) -> PyResult<Quaternion> {
        Ok(Quaternion {
            inner: self.inner.conjugate(),
        })
    }

    fn __mul__(&self, other: &PyAny) -> PyResult<PyObject> {
        // Multiply quaternion by quaternion
        if other.is_instance_of::<Quaternion>().unwrap() {
            let q: PyRef<Quaternion> = other.extract()?;
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                return Ok(Quaternion {
                    inner: self.inner * q.inner,
                }
                .into_py(py));
            })
        }
        // This incorrectly matches for all PyArray types
        else if other.is_instance_of::<np::PyArray2<f64>>().unwrap() {
            // So, check for 2D condition
            match other.extract::<np::PyReadonlyArray2<f64>>() {
                Ok(v) => {
                    if v.dims()[1] != 3 {
                        return Err(pyo3::exceptions::PyTypeError::new_err(
                            "Invalid rhs. 2nd dimension must be 3 in size",
                        ));
                    }
                    let rot = self.inner.to_rotation_matrix();
                    let qmat = rot.matrix().conjugate();

                    pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                        unsafe {
                            let nd = np::ndarray::ArrayView2::from_shape_ptr((3, 3), qmat.as_ptr());
                            let res = v.as_array().dot(&nd).to_pyarray(py);

                            Ok(res.into_py(py))
                        }
                    })
                }
                // If not, check for 1D condition
                Err(_) => match other.extract::<np::PyReadonlyArray1<f64>>() {
                    Ok(v1) => {
                        if v1.len() != 3 {
                            return Err(pyo3::exceptions::PyValueError::new_err(
                                "rhs 1D array must have 3 elements",
                            ));
                        }
                        let mut s = Vec3::zeros();

                        pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                            unsafe {
                                std::ptr::copy_nonoverlapping(
                                    v1.as_raw_array().as_ptr(),
                                    s.as_mut_ptr(),
                                    3,
                                );

                                let vout = self.inner * s;
                                let vnd = np::PyArray1::<f64>::from_vec(
                                    py,
                                    vec![vout[0], vout[1], vout[2]],
                                );
                                Ok(vnd.into_py(py))
                            }
                        })
                    }
                    // Input is incorrect size...
                    Err(_) => {
                        return Err(pyo3::exceptions::PyIndexError::new_err(
                            "RHS must be 1x3 or nx3",
                        ));
                    }
                },
            }
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err("Invalid rhs"))
        }
    }
}
