use nalgebra as na;
use numpy as np;
use numpy::ndarray::array;
use pyo3::prelude::*;
use pyo3::{exceptions, Python};

type Quat = na::UnitQuaternion<f64>;
type Vec3 = na::Vector3<f64>;

#[pyclass]
pub struct Quaternion {
    inner: Quat,
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
    fn from_axis_angle(axis: np::PyReadonlyArray1<f64>, angle: f64) -> PyResult<Self> {
        let v: Vec3 = Vec3::new(axis.as_array()[0], axis.as_array()[1], axis.as_array()[2]);
        let u = na::UnitVector3::try_new(v, 1.0e-9);
        match u {
            Some(ax) => Ok(Quaternion {
                inner: Quat::from_axis_angle(&ax, angle),
            }),
            None => {
                let err = exceptions::PyArithmeticError::new_err("Axis norm is 0");
                Err(err)
            }
        }
    }

    #[getter]
    fn angle(&self) -> PyResult<f64> {
        Ok(self.inner.angle())
    }

    #[getter]
    fn axis(&self) -> PyResult<Py<numpy::PyArray1<f64>>> {
        let a = match self.inner.axis() {
            Some(ax) => ax,
            None => Vec3::x_axis(),
        };
        let gil = Python::acquire_gil();
        let pyarray = np::PyArray1::from_array(gil.python(), &array![a[0], a[1], a[2]]);

        Ok(pyarray.to_owned())
    }
}
