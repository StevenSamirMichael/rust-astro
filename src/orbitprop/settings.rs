//! Orbit Propagation Settings

pub struct PropSettings {
    pub gravity_order: u16,
    pub gravity_interp_dt_secs: f64,
    pub moon_interp_dt_secs: f64,
    pub sun_interp_dt_secs: f64,
    pub abs_error: f64,
    pub rel_error: f64,
    pub dt_secs: f64,
}

impl PropSettings {
    pub fn new() -> PropSettings {
        PropSettings {
            gravity_order: 6,
            gravity_interp_dt_secs: 120.0,
            moon_interp_dt_secs: 120.0,
            sun_interp_dt_secs: 120.0,
            abs_error: 1e-8,
            rel_error: 1e-8,
            dt_secs: 1.0,
        }
    }
}
