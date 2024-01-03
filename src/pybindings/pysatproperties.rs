use crate::orbitprop::SatPropertiesStatic;

use pyo3::prelude::*;

#[pyclass(name = "satproperties_static")]
#[derive(Clone, Debug)]
pub struct PySatProperties {
    pub inner: SatPropertiesStatic,
}

#[pymethods]
impl PySatProperties {
    #[new]
    fn py_new() -> PyResult<Self> {
        Ok(PySatProperties {
            inner: SatPropertiesStatic::new(0.0, 0.0),
        })
    }

    #[getter]
    fn get_craoverm(&self) -> f64 {
        self.inner.craoverm
    }

    #[getter]
    fn get_cdaoverm(&self) -> f64 {
        self.inner.cdaoverm
    }

    #[setter]
    fn set_craoverm(&mut self, craoverm: f64) -> PyResult<()> {
        self.inner.craoverm = craoverm;
        Ok(())
    }

    #[setter]
    fn set_cdaoverm(&mut self, cdaoverm: f64) -> PyResult<()> {
        self.inner.cdaoverm = cdaoverm;
        Ok(())
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner.to_string())
    }
}
