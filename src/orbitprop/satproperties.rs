use super::propagator::SimpleState;
use crate::AstroTime;

pub trait SatProperties {
    // Coefficient of drag times normal area over mass
    fn cd_a_over_m(&self, tm: &AstroTime, state: &SimpleState) -> f64;

    // Coefficient of radiation pressure times normal area over mass
    fn cr_a_over_m(&self, tm: &AstroTime, state: &SimpleState) -> f64;
}
