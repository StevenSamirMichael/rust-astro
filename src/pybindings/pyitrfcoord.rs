use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::types::PyList;

use numpy::{PyArray1, PyReadonlyArray1};

use crate::itrfcoord::ITRFCoord;
use crate::types::Vec3;

use super::Quaternion;

///
/// Representation of a coordinate in the
/// International Terrestrial Reference Frame (ITRF)
///
/// This coordinate object can be created from and also
/// output to Geodetic coordinates (latitude, longitude,
/// height above ellipsoid).
///
/// Functions are also available to provide rotation
/// quaternions to the East-North-Up frame
/// and North-East-Down frame at this coordinate
///
#[pyclass(name = "itrfcoord")]
pub struct PyITRFCoord {
    pub inner: ITRFCoord,
}

use std::f64::consts::PI;
const DEG2RAD: f64 = PI / 180.;

#[pymethods]
impl PyITRFCoord {
    // Represent a coordinate in the ITRF (International Terrestrial Reference Frame)
    #[new]
    #[pyo3(signature=(x = 0.0, y = 0.0, z = 0.0))]
    fn new(x: f64, y: f64, z: f64) -> PyResult<Self> {
        Ok(PyITRFCoord {
            inner: ITRFCoord::from_slice(&[x, y, z]).unwrap(),
        })
    }

    // From a 3-element vector or list
    #[staticmethod]
    fn from_vector(v: &PyAny) -> PyResult<Self> {
        if v.is_instance_of::<PyList>() {
            let lv = v.extract::<Vec<f64>>()?;
            if lv.len() != 3 {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "Invalid number of elements",
                ));
            }
            Ok(PyITRFCoord {
                inner: ITRFCoord::from_slice(lv.as_slice()).unwrap(),
            })
        } else if v.is_instance_of::<PyArray1<f64>>() {
            let nv = v.extract::<PyReadonlyArray1<f64>>()?;
            if nv.len() != 3 {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "Invalid number of elements",
                ));
            }
            Ok(PyITRFCoord {
                inner: ITRFCoord::from_slice(nv.as_array().as_slice().unwrap()).unwrap(),
            })
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err("Invalid input type"));
        }
    }

    ///
    /// Create coordinate from input geodetic
    /// Optional inputs, in order:
    ///
    /// latitude, radians
    /// longitude, radians
    /// heigth above ellipsoid, meters
    ///
    ///
    /// Optional kwargs:
    ///
    /// latitude_deg: latitude, degrees
    /// longitude_deg: longitude, degrees
    /// latitude_rad: latitude, radians
    /// longitude_rad: longitude, radians
    /// altitude: height above ellipsoid, meters
    ///
    #[staticmethod]
    #[pyo3(signature=(latitude=0.0, longitude=0.0, altitude=0.0, **kwargs))]
    fn from_geodetic(
        mut latitude: f64,
        mut longitude: f64,
        mut altitude: f64,
        kwargs: Option<&PyDict>,
    ) -> PyResult<Self> {
        if kwargs.is_some() {
            let kw = kwargs.unwrap();

            match kw.get_item("latitude_deg")? {
                Some(v) => {
                    latitude = v.extract::<f64>()? * DEG2RAD;
                }
                None => (),
            }
            match kw.get_item("longitude_deg")? {
                Some(v) => {
                    longitude = v.extract::<f64>()? * DEG2RAD;
                }
                None => (),
            }

            match kw.get_item("latitude_rad")? {
                Some(v) => {
                    latitude = v.extract::<f64>()?;
                }
                None => (),
            }
            match kw.get_item("longitude_rad")? {
                Some(v) => {
                    longitude = v.extract::<f64>()?;
                }
                None => (),
            }
            match kw.get_item("altitude")? {
                Some(v) => {
                    altitude = v.extract::<f64>()?;
                }
                None => (),
            }
            match kw.get_item("height")? {
                Some(v) => {
                    altitude = v.extract::<f64>()?;
                }
                None => (),
            }
        }

        Ok(PyITRFCoord {
            inner: ITRFCoord::from_geodetic_rad(latitude, longitude, altitude),
        })
    }

    #[getter]
    /// Latitude in degrees
    fn get_latitude_deg(&self) -> f64 {
        self.inner.latitude_deg()
    }

    /// Longitude in degrees
    #[getter]
    fn get_longitude_deg(&self) -> f64 {
        self.inner.longitude_deg()
    }

    /// Latitude in radians
    #[getter]
    fn get_latitude_rad(&self) -> f64 {
        self.inner.latitude_rad()
    }

    /// Longitude in radians
    #[getter]
    fn get_longitude_rad(&self) -> f64 {
        self.inner.longitude_rad()
    }

    /// Height above ellipsoid in meters
    #[getter]
    fn get_height(&self) -> f64 {
        self.inner.hae()
    }

    /// Height above ellipsoid, meters
    #[getter]
    fn get_altitude(&self) -> f64 {
        self.inner.hae()
    }

    /// Return Tuple with latitude in rad, longitude in rad,
    /// height above ellipsoid in meters
    #[getter]
    fn get_geodetic_rad(&self) -> (f64, f64, f64) {
        self.inner.to_geodetic_rad()
    }

    /// Return tuple with latitude in deg, longitude in deg,
    /// height above ellipsoid in meters
    #[getter]
    fn get_geodetic_deg(&self) -> (f64, f64, f64) {
        self.inner.to_geodetic_deg()
    }

    /// Return vector representing ITRF Cartesian coordinate
    /// in meters
    #[getter]
    fn get_vec(&self) -> PyObject {
        pyo3::Python::with_gil(|py| -> PyObject {
            numpy::PyArray::from_slice(py, self.inner.itrf.data.as_slice()).to_object(py)
        })
    }

    fn __str__(&self) -> String {
        let (lat, lon, hae) = self.inner.to_geodetic_deg();
        format!(
            "ITRFCoord(lat: {:8.4} deg, lon: {:8.4} deg, hae: {:5.2} km)",
            lat,
            lon,
            hae / 1.0e3
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }

    /// Quaternion representing rotation from
    /// North-East-Down (NED) coordinate frame
    /// to International Terrestrial Reference Frame
    /// (ITRF) at this point
    fn qned2itrf(&self) -> Quaternion {
        Quaternion {
            inner: self.inner.q_ned2itrf(),
        }
    }

    /// Quaternion representing rotation from
    /// East-North-Up (ENU) coordinate frame
    /// to International Terrestrial Reference Frame
    /// (ITRF) at this point
    fn qenu2itrf(&self) -> Quaternion {
        Quaternion {
            inner: self.inner.q_ned2itrf(),
        }
    }

    /// Return East-North-Up location of input
    /// coordinate relative to self
    fn to_enu(&self, other: &Self) -> PyObject {
        let v: Vec3 = self.inner.q_enu2itrf().conjugate() * (self.inner.itrf - other.inner.itrf);
        pyo3::Python::with_gil(|py| -> PyObject {
            numpy::PyArray::from_slice(py, v.data.as_slice()).to_object(py)
        })
    }

    /// Return North-East-Down location of input
    /// coordinate relative to self
    fn to_ned(&self, other: &Self) -> PyObject {
        let v: Vec3 = self.inner.q_ned2itrf().conjugate() * (self.inner.itrf - other.inner.itrf);
        pyo3::Python::with_gil(|py| -> PyObject {
            numpy::PyArray::from_slice(py, v.data.as_slice()).to_object(py)
        })
    }
}

impl IntoPy<PyObject> for ITRFCoord {
    fn into_py(self, py: Python<'_>) -> PyObject {
        PyITRFCoord { inner: self }.into_py(py)
    }
}

impl<'b> From<&'b PyITRFCoord> for &'b ITRFCoord {
    fn from<'a>(s: &'a PyITRFCoord) -> &'a ITRFCoord {
        &s.inner
    }
}
