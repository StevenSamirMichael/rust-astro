use cpython::{self, py_class, PyObject, PyResult, PyString, ToPyObject};
use nalgebra as na;
type Quat = na::UnitQuaternion<f64>;
type Vec3 = na::Vector3<f64>;
use std::cell::RefCell;

py_class! {
    pub class Quaternion |py| {
        data quat: RefCell<Quat>;

        def __new__(_cls) -> PyResult<Quaternion>
        {
            Self::create_instance(py, RefCell::new({
                    Quat::new(Vec3::z() * 0.0)
                }))
        }

        def __str__(&self) -> PyResult<PyString>
        {
            let q = self.quat(py).borrow();
            let a = match q.axis() {
                Some(ax) => ax,
                None => Vec3::x_axis()
            };
            let f = format!("Quaternion(Axis = [{:7.4}, {:7.4}, {:7.4}] , Angle = {:7.4} rad)", a[0], a[1], a[2], q.angle());
            Ok(PyString::new(py, f.as_str()))
        }

        def conjugate(&self) -> PyResult<Quaternion>
        {
            Self::create_instance(py, RefCell::new(self.quat(py).borrow().conjugate()))
        }

        def conj(&self) -> PyResult<Quaternion>
        {
            Self::create_instance(py, RefCell::new(self.quat(py).borrow().conjugate()))
        }


        def __mul__(lhs, rhs) -> PyResult<impl ToPyObject>
        {

            let h: &PyObject = lhs;
            let t = h.get_type(py);
            let q: Quaternion = h.extract(py)?;


            t.is_instance(py, h);


            let c = Quat::from_axis_angle(&Vec3::x_axis(), 1.0);
            Self::create_instance(py, RefCell::new(c))
        }


        @staticmethod def rotx(theta: f64)->PyResult<Quaternion>
        {
            let c = Quat::from_axis_angle(&Vec3::x_axis(), theta);
            Self::create_instance(py, RefCell::new(c))
        }

    }
}
