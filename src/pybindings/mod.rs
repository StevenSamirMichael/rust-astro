use cpython::py_module_initializer;

mod pyquaternion;

py_module_initializer! {
    astro, |py, module| {
        module.add(py, "__doc__", "This module is implemented in Rust.")?;

        module.add_class::<pyquaternion::Quaternion>(py)?;

        Ok(())
    }
}
