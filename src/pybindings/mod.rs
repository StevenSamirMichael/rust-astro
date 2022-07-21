use pyo3::prelude::*;

mod pyquaternion;

use pyquaternion::Quaternion;

#[pyfunction]
fn sum(a: i32, b: i32) -> PyResult<i32> {
    Ok(a + b)
}

#[pymodule]
pub fn astro(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(sum))?;

    m.add_class::<Quaternion>()?;
    Ok(())
}
