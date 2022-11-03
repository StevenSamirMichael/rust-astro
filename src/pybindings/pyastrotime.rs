use pyo3::prelude::*;
use pyo3::types::PyDateTime;

use crate::astrotime::{self, AstroTime, Scale, TimeInput, TimeInputType};

use numpy as np;
use std::cell::RefCell;

#[pyclass(name = "timescale")]
pub enum PyTimeScale {
    /// Invalid time scale
    Invalid = Scale::INVALID as isize,
    /// Universal Time Coordinate
    UTC = Scale::UTC as isize,
    /// Terrestrial Time
    TT = Scale::TT as isize,
    /// UT1
    UT1 = Scale::UT1 as isize,
    /// International Atomic Time
    TAI = Scale::TAI as isize,
    /// Global Positioning System (GPS) Time
    GPS = Scale::GPS as isize,
    /// Barycentric Dynamical Time
    TDB = Scale::TDB as isize,
}

impl From<&PyTimeScale> for astrotime::Scale {
    fn from(s: &PyTimeScale) -> astrotime::Scale {
        match s {
            PyTimeScale::Invalid => Scale::INVALID,
            PyTimeScale::UTC => Scale::UTC,
            PyTimeScale::TT => Scale::TT,
            PyTimeScale::UT1 => Scale::UT1,
            PyTimeScale::TAI => Scale::TAI,
            PyTimeScale::GPS => Scale::GPS,
            PyTimeScale::TDB => Scale::TDB,
        }
    }
}

impl IntoPy<PyObject> for astrotime::Scale {
    fn into_py(self, py: Python<'_>) -> PyObject {
        let ts: PyTimeScale = match self {
            Scale::INVALID => PyTimeScale::Invalid,
            Scale::UTC => PyTimeScale::UTC,
            Scale::TT => PyTimeScale::TT,
            Scale::UT1 => PyTimeScale::UT1,
            Scale::TAI => PyTimeScale::TAI,
            Scale::GPS => PyTimeScale::GPS,
            Scale::TDB => PyTimeScale::TDB,
        };
        ts.into_py(py)
    }
}

#[pyclass(name = "time")]
#[derive(PartialEq, PartialOrd, Copy, Clone, Debug)]
pub struct PyAstroTime {
    pub inner: AstroTime,
}

#[pymethods]
impl PyAstroTime {
    #[new]
    fn py_new() -> PyResult<Self> {
        match AstroTime::now() {
            Ok(v) => Ok(PyAstroTime { inner: v }),
            Err(_) => Err(pyo3::exceptions::PyOSError::new_err(
                "Could not get current time",
            )),
        }
    }

    /// Return current time
    #[staticmethod]
    fn now() -> PyResult<Self> {
        match AstroTime::now() {
            Ok(v) => Ok(PyAstroTime { inner: v }),
            Err(_) => Err(pyo3::exceptions::PyOSError::new_err(
                "Could not get current time",
            )),
        }
    }

    /// Return time object representing input
    /// Gregorian year, month (1=January, 2=February, ...), and
    /// day of month, beginning with 1.  Inputs assumed to be UTC
    #[staticmethod]
    fn from_date(year: u32, month: u32, day: u32) -> PyResult<Self> {
        Ok(PyAstroTime {
            inner: AstroTime::from_date(year, month, day),
        })
    }

    /// Convert time object to UTC Gegorian date, with
    /// returns tuple with 3 elements:
    /// 1 : Gregorian Year
    /// 2 : Gregorian month (1 = January, 2 = February, ...)
    /// 3 : Day of month, beginning with 1
    ///
    fn to_date(&self) -> (u32, u32, u32) {
        self.inner.to_date()
    }

    /// Convert time object to UTC Gegorian date and time, with
    /// returns tuple with 6 elements:
    /// 1 : Gregorian Year
    /// 2 : Gregorian month (1 = January, 2 = February, ...)
    /// 3 : Day of month, beginning with 1
    /// 4 : Hour of day, in range [0,23]
    /// 5 : Minute of hour, in range [0,59]
    /// 6 : floating point second of minute, in range [0,60)
    ///
    fn to_gregorian(&self) -> (u32, u32, u32, u32, u32, f64) {
        self.inner.to_datetime()
    }

    /// Convert UTC Gegorian date and time to time object with
    /// 6-element input:
    /// 1 : Gregorian Year
    /// 2 : Gregorian month (1 = January, 2 = February, ...)
    /// 3 : Day of month, beginning with 1
    /// 4 : Hour of day, in range [0,23]
    /// 5 : Minute of hour, in range [0,59]
    /// 6 : floating point second of minute, in range [0,60)
    ///
    #[staticmethod]
    fn from_gregorian(
        year: u32,
        month: u32,
        day: u32,
        hour: u32,
        min: u32,
        sec: f64,
    ) -> PyResult<Self> {
        Ok(PyAstroTime {
            inner: AstroTime::from_datetime(year, month, day, hour, min, sec),
        })
    }

    /// Convert to Python datetime object
    fn datetime(&self) -> PyResult<PyObject> {
        pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
            let timestamp: f64 = self.to_unixtime();
            Ok(PyDateTime::from_timestamp(py, timestamp, None)?.into_py(py))
        })
    }

    fn to_mjd(&self, scale: &PyTimeScale) -> f64 {
        self.inner.to_mjd(scale.into())
    }

    fn to_jd(&self, scale: &PyTimeScale) -> f64 {
        self.inner.to_jd(scale.into())
    }

    fn to_unixtime(&self) -> f64 {
        self.inner.to_unixtime()
    }

    fn __add__(&self, other: &PyAny) -> PyResult<PyObject> {
        // Numpy array of floats
        if other.is_instance_of::<np::PyArray1<f64>>().unwrap() {
            let parr: np::PyReadonlyArray1<f64> = other.extract().unwrap();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                let objarr = parr
                    .as_array()
                    .into_iter()
                    .map(|x| {
                        let obj = PyAstroTime {
                            inner: self.inner + *x,
                        };
                        obj.into_py(py)
                    })
                    .into_iter();
                let parr = np::PyArray1::<PyObject>::from_iter(py, objarr);
                Ok(parr.into_py(py))
            })
        }
        // list of floats
        else if other.is_instance_of::<pyo3::types::PyList>().unwrap() {
            let v = other.extract::<Vec<f64>>().unwrap();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                let objarr = v
                    .into_iter()
                    .map(|x| {
                        let pyobj = PyAstroTime {
                            inner: self.inner + x,
                        };
                        pyobj.into_py(py)
                    })
                    .into_iter();

                let parr = np::PyArray1::<PyObject>::from_iter(py, objarr);
                Ok(parr.into_py(py))
            })
        }
        // Constant number
        else if other.is_instance_of::<pyo3::types::PyFloat>().unwrap()
            || other.is_instance_of::<pyo3::types::PyInt>().unwrap()
            || other.is_instance_of::<pyo3::types::PyLong>().unwrap()
        {
            let dt: f64 = other.extract::<f64>().unwrap();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                Ok(PyAstroTime {
                    inner: self.inner + dt,
                }
                .into_py(py))
            })
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err(
                "Invalid type for rhs",
            ))
        }
    }

    fn __sub__(&self, other: &PyAny) -> PyResult<PyObject> {
        // Numpy array of floats
        if other.is_instance_of::<np::PyArray1<f64>>().unwrap() {
            let parr: np::PyReadonlyArray1<f64> = other.extract().unwrap();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                let objarr = parr
                    .as_array()
                    .into_iter()
                    .map(|x| {
                        let obj = PyAstroTime {
                            inner: self.inner - *x,
                        };
                        obj.into_py(py)
                    })
                    .into_iter();
                let parr = np::PyArray1::<PyObject>::from_iter(py, objarr);
                Ok(parr.into_py(py))
            })
        }
        // list of floats
        else if other.is_instance_of::<pyo3::types::PyList>().unwrap() {
            let v = other.extract::<Vec<f64>>().unwrap();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                let objarr = v
                    .into_iter()
                    .map(|x| {
                        let pyobj = PyAstroTime {
                            inner: self.inner - x,
                        };
                        pyobj.into_py(py)
                    })
                    .into_iter();

                let parr = np::PyArray1::<PyObject>::from_iter(py, objarr);
                Ok(parr.into_py(py))
            })
        }
        // Constant number
        else if other.is_instance_of::<pyo3::types::PyFloat>().unwrap()
            || other.is_instance_of::<pyo3::types::PyInt>().unwrap()
            || other.is_instance_of::<pyo3::types::PyLong>().unwrap()
        {
            let dt: f64 = other.extract::<f64>().unwrap();
            pyo3::Python::with_gil(|py| -> PyResult<PyObject> {
                Ok(PyAstroTime {
                    inner: self.inner - dt,
                }
                .into_py(py))
            })
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err(
                "Invalid type for rhs",
            ))
        }
    }

    fn add_utc_days(&self, days: f64) -> PyAstroTime {
        PyAstroTime {
            inner: self.inner.add_utc_days(days),
        }
    }

    fn __str__(&self) -> PyResult<String> {
        let dt = self.inner.to_datetime();
        Ok(format!(
            "{:04}-{:02}-{:02} {:02}:{:02}:{:06.3}Z",
            dt.0, dt.1, dt.2, dt.3, dt.4, dt.5
        ))
    }

    fn __repr__(&self) -> PyResult<String> {
        self.__str__()
    }
}

impl IntoPy<PyObject> for astrotime::AstroTime {
    fn into_py(self, py: Python<'_>) -> PyObject {
        let ts: PyAstroTime = PyAstroTime { inner: self };
        ts.into_py(py)
    }
}

impl<'b> From<&'b PyAstroTime> for &'b astrotime::AstroTime {
    fn from<'a>(s: &'a PyAstroTime) -> &'a astrotime::AstroTime {
        &s.inner
    }
}

impl TimeInputType for &PyAny {
    fn to_time_input(&self) -> TimeInput {
        // "Scalar" time input case
        if self.is_instance_of::<PyAstroTime>().unwrap() {
            let tm: PyAstroTime = self.extract().unwrap();
            TimeInput::Single(tm.inner)
        }
        // List case
        else if self.is_instance_of::<pyo3::types::PyList>().unwrap() {
            match self.extract::<Vec<PyAstroTime>>() {
                Ok(v) => TimeInput::Array(RefCell::new(v.into_iter().map(|x| x.inner).collect())),
                Err(e) => TimeInput::Error(e.to_string()),
            }
        }
        // numpy array case
        else if self.is_instance_of::<numpy::PyArray1<PyObject>>().unwrap() {
            match self.extract::<numpy::PyReadonlyArray1<PyObject>>() {
                Ok(v) => pyo3::Python::with_gil(|py| -> TimeInput {
                    // Extract times from numpya array of objects
                    let tmarray: Result<Vec<AstroTime>, _> = v
                        .as_array()
                        .into_iter()
                        .map(|p| -> Result<AstroTime, _> {
                            match p.extract::<PyAstroTime>(py) {
                                Ok(v2) => Ok(v2.inner),
                                Err(v2) => Err(v2),
                            }
                        })
                        .collect();
                    if !tmarray.is_ok() {
                        return TimeInput::Error(String::from("Invalid data type in numpy array"));
                    }
                    TimeInput::Array(RefCell::new(tmarray.unwrap()))
                }),

                Err(e) => TimeInput::Error(e.to_string()),
            }
        } else {
            TimeInput::Error(String::from(
                "Input is not astro::time or ndarray or list of astro::time",
            ))
        }
    }
}
