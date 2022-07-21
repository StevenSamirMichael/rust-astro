use nalgebra as na;
use nalgebra::{Dynamic, Matrix, SliceStorageMut};
use numpy as np;
use numpy::ndarray::array;
use numpy::npyffi::objects::PyArrayObject;
use numpy::ToPyArray;
use pyo3::prelude::*;
//use pyo3::{exceptions, Python};

type Quat = na::UnitQuaternion<f64>;
type Vec3 = na::Vector3<f64>;

#[pyclass]
#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Quaternion {
    inner: Quat,
}

fn matrix_to_numpy<'py, N, R, C, S>(
    py: pyo3::Python<'py>,
    matrix: &na::Matrix<N, R, C, S>,
) -> pyo3::PyObject
where
    N: nalgebra::Scalar + numpy::Element,
    R: nalgebra::Dim,
    C: nalgebra::Dim,
    S: nalgebra::storage::Storage<N, R, C>,
{
    unsafe {
        let array = np::PyArray::new(py, (matrix.nrows(), matrix.ncols()), false);

        for r in 0..matrix.nrows() {
            for c in 0..matrix.ncols() {
                *array.uget_mut((r, c)) = matrix[(r, c)].clone();
            }
        }
        array.into_py(py)
    }
}

#[allow(clippy::missing_safety_doc)]
unsafe fn matrix_slice_mut_from_numpy_ptr<'a, N, R, C>(
    array: *mut PyArrayObject,
) -> nalgebra::MatrixSliceMut<'a, N, R, C, Dynamic, Dynamic>
where
    N: nalgebra::Scalar + numpy::Element,
    R: nalgebra::Dim,
    C: nalgebra::Dim,
{
    //let array = cast_to_py_array(array)?;

    let input_rows = *(*array).dimensions.add(0) as usize;
    let input_cols = *(*array).dimensions.add(1) as usize;
    let shape = (R::from_usize(input_rows), C::from_usize(input_cols));
    //check_array_alignment(array)?;

    let row_stride = Dynamic::new(*(*array).strides.add(0) as usize / std::mem::size_of::<N>());
    let col_stride = Dynamic::new(*(*array).strides.add(1) as usize / std::mem::size_of::<N>());
    let storage = SliceStorageMut::<N, R, C, Dynamic, Dynamic>::from_raw_parts(
        (*array).data as *mut N,
        shape,
        (row_stride, col_stride),
    );

    Matrix::from_data(storage)
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
        Ok(array![a[0], a[1], a[2]].to_pyarray(gil.python()).to_owned())
    }

    fn __mul__(&self, other: &PyAny) -> PyResult<PyObject> {
        if other.is_instance_of::<Quaternion>().unwrap() {
            let q: PyRef<Quaternion> = other.extract()?;
            let gil = Python::acquire_gil();
            return Ok(Quaternion {
                inner: self.inner * q.inner,
            }
            .into_py(gil.python()));
        } else if other.is_instance_of::<np::PyArray1<f64>>().unwrap() {
            let v: &np::PyArray1<f64> = other.extract()?;
            if v.len() != 3 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "rhs 1D array must have 3 elements",
                ));
            }
            let s = Vec3::new(
                v.get_owned(0).unwrap(),
                v.get_owned(1).unwrap(),
                v.get_owned(2).unwrap(),
            );
            let vout = self.inner * s;
            let gil = Python::acquire_gil();
            let vnd = np::PyArray1::<f64>::from_vec(gil.python(), vec![vout[0], vout[1], vout[2]]);
            return Ok(vnd.into_py(gil.python()));
        } else if other.is_instance_of::<np::PyArray2<f64>>().unwrap() {
            let v: &np::PyArray2<f64> = other.extract()?;
            if v.dims()[1] != 3 {
                return Err(pyo3::exceptions::PyTypeError::new_err("Invalid rhs"));
            }
            let qmat = self.inner.to_rotation_matrix();

            Err(pyo3::exceptions::PyTypeError::new_err(
                "Not implemented yet",
            ))
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err("Invalid rhs"))
        }
    }
}
