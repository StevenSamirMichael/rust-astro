pub struct RKAdaptiveSettings {
    pub abserror: f64,
    pub relerror: f64,
    pub dense_output: bool,
}

impl RKAdaptiveSettings {
    pub fn default() -> RKAdaptiveSettings {
        RKAdaptiveSettings {
            abserror: 1.0e-10,
            relerror: 1.0e-6,
            dense_output: false,
        }
    }
}
