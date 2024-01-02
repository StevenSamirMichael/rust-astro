use super::pyastrotime::PyAstroTime;
use super::pypropsettings::PyPropSettings;

use nalgebra as na;
use numpy as np;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyString};

fn kwargs_or_default<'a, T>(kwargs: &mut Option<&'a PyDict>, name: &str, default: T) -> PyResult<T>
where
    T: FromPyObject<'a>,
{
    match kwargs.is_some() {
        true => {
            let kw = kwargs.unwrap();
            match kw.get_item(name)? {
                None => Ok(default),
                Some(v) => {
                    kw.del_item(name)?;
                    Ok(v.extract::<T>()?)
                }
            }
        }
        false => Ok(default),
    }
}

fn kwargs_or_none<'a, T>(kwargs: &mut Option<&'a PyDict>, name: &str) -> PyResult<Option<T>>
where
    T: FromPyObject<'a>,
{
    match kwargs.is_some() {
        true => {
            let kw = kwargs.unwrap();
            match kw.get_item(name)? {
                None => Ok(None),
                Some(v) => {
                    kw.del_item(name)?;
                    Ok(Some(v.extract::<T>()?))
                }
            }
        }
        false => Ok(None),
    }
}

#[pyfunction(signature=(pos, vel, start, **kwargs))]
pub fn propagate(
    pos: &np::PyArray1<f64>,
    vel: &np::PyArray1<f64>,
    start: &PyAstroTime,
    kwargs: Option<&PyDict>,
) -> PyResult<Py<PyAny>> {
    let mut mkwargs = kwargs.clone();

    if pos.len() != 3 || vel.len() != 3 {
        return Err(pyo3::exceptions::PyRuntimeError::new_err(
            "Position and velocity must be 1-d numpy arrays with length 3",
        ));
    }
    let pypropsettings: Option<PyPropSettings> = kwargs_or_none(&mut mkwargs, "propsettings")?;

    let propsettings = match pypropsettings {
        Some(p) => p.inner,
        None => crate::orbitprop::PropSettings::default(),
    };
    let mut dt_secs: Option<f64> = kwargs_or_none(&mut mkwargs, "dt_secs")?;
    let dt_days: Option<f64> = kwargs_or_none(&mut mkwargs, "dt_days")?;
    let duration_secs: Option<f64> = kwargs_or_none(&mut mkwargs, "duration_secs")?;
    let duration_days: Option<f64> = kwargs_or_none(&mut mkwargs, "duration_days")?;
    let pystoptime: Option<PyAstroTime> = kwargs_or_none(&mut mkwargs, "stoptime")?;
    //let output_phi: bool = kwargs_or_default(&kwargs, "output_phi", false)?;

    match dt_days {
        None => (),
        Some(v) => dt_secs = Some(v * 86400.0),
    };

    // Look for extraneous kwargs and return error
    if mkwargs.is_some() {
        if !mkwargs.unwrap().is_empty() {
            let keystring: String =
                mkwargs
                    .unwrap()
                    .iter()
                    .fold(String::from(""), |acc, (k, _v)| {
                        let mut a2 = acc.clone();
                        a2.push_str(k.downcast::<PyString>().unwrap().to_str().unwrap());
                        a2.push_str(", ");
                        a2
                    });
            let s = format!("Invalid kwargs: {}", keystring);
            return Err(pyo3::exceptions::PyRuntimeError::new_err(s));
        }
    }

    if duration_days == None && pystoptime == None && duration_secs == None {
        return Err(pyo3::exceptions::PyRuntimeError::new_err(
            "Must set either duration or stop time",
        ));
    }
    let stoptime = match pystoptime {
        Some(p) => p.inner,
        None => {
            start.inner
                + match duration_days {
                    Some(v) => v,
                    None => duration_secs.unwrap() / 86400.0,
                }
        }
    };

    // Create the state to propagate
    let mut pv = na::Vector6::<f64>::zeros();
    pv.fixed_view_mut::<3, 1>(0, 0)
        .copy_from_slice(unsafe { pos.as_slice().unwrap() });
    pv.fixed_view_mut::<3, 1>(3, 0)
        .copy_from_slice(unsafe { vel.as_slice().unwrap() });

    // Finally, do the propagation
    let res = match crate::orbitprop::propagate(
        &pv,
        &start.inner,
        &stoptime,
        dt_secs,
        &propsettings,
        None,
    ) {
        Ok(v) => v,
        Err(e) => {
            let estring = format!("Error propagating: {}", e.to_string());
            return Err(pyo3::exceptions::PyRuntimeError::new_err(estring));
        }
    };

    pyo3::Python::with_gil(|py| -> PyResult<Py<PyAny>> {
        let r = PyDict::new(py);

        let d = PyDict::new(py);
        d.set_item("num_eval", res.num_eval)?;
        d.set_item("accepted_steps", res.accepted_steps)?;
        d.set_item("rejected_steps", res.rejected_steps)?;

        let tm: Vec<Py<PyAny>> = res
            .time
            .iter()
            .map(|x| PyAstroTime { inner: x.clone() }.into_py(py))
            .collect();

        let n = res.state.len();
        let pos = unsafe { np::PyArray2::<f64>::new(py, [n, 3], false) };
        for idx in 0..n {
            unsafe {
                std::ptr::copy_nonoverlapping(
                    res.state[idx].as_ptr(),
                    pos.as_raw_array_mut().as_mut_ptr().offset(idx as isize * 3),
                    3,
                );
            }
        }
        let vel = unsafe { np::PyArray2::<f64>::new(py, [n, 3], false) };
        for idx in 0..n {
            unsafe {
                std::ptr::copy_nonoverlapping(
                    res.state[idx].as_ptr().offset(3),
                    vel.as_raw_array_mut().as_mut_ptr().offset(idx as isize * 3),
                    3,
                );
            }
        }
        r.set_item("stats", d)?;
        r.set_item("time", tm)?;
        r.set_item("pos", pos)?;
        r.set_item("vel", vel)?;

        Ok(r.to_object(py))
    })
}
