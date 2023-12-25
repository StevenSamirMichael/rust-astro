//! Orbit Propagation Settings

pub struct PropSettings {
    pub gravity_order: u16,
    pub gravity_interp_dt_secs: f64,
    pub moon_interp_dt_secs: f64,
    pub sun_interp_dt_secs: f64,
    pub abs_error: f64,
    pub rel_error: f64,
    pub use_spaceweather: bool,
}

impl PropSettings {
    pub fn new() -> PropSettings {
        PropSettings {
            gravity_order: 6,
            gravity_interp_dt_secs: 60.0,
            moon_interp_dt_secs: 60.0,
            sun_interp_dt_secs: 60.0,
            abs_error: 1e-10,
            rel_error: 1e-10,
            use_spaceweather: true,
        }
    }
}
