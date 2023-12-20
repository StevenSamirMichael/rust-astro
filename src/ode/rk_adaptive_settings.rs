pub struct RKAdaptiveSettings {
    pub abserror: f64,
    pub relerror: f64,
    pub minfac: f64,
    pub maxfac: f64,
    pub safetyfac: f64,
    pub gamma: f64,
    pub dtmin: f64,
}

impl RKAdaptiveSettings {
    pub fn default() -> RKAdaptiveSettings {
        RKAdaptiveSettings {
            abserror: 1.0e-8,
            relerror: 1.0e-8,
            minfac: 0.2,
            maxfac: 5.0,
            safetyfac: 0.9,
            gamma: 0.9,
            dtmin: 1.0e-6,
        }
    }
}
