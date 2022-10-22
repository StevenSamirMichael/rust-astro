use super::force_gravity::ForceGravity;
use super::force_moon::ForceMoon;
use super::force_sun::ForceSun;

use super::forceterm::*;
use super::settings::PropSettings;
use super::types::*;

use crate::astrotime::AstroTime;
use crate::utils::{astroerr, AstroResult};
use ode_solvers::dopri5::*;
use std::vec::Vec;

#[derive(Debug)]
pub struct PropagationResult<T> {
    pub time: Vec<AstroTime>,
    pub state: Vec<T>,
    pub accepted_steps: u32,
    pub rejected_steps: u32,
    pub num_eval: u32,
}

struct Propagation<'a> {
    start: AstroTime,
    forces: Vec<&'a dyn ForceTerm<SimpleState>>,
}

impl<'a> ode_solvers::System<SimpleState> for Propagation<'a> {
    fn system(&self, x: f64, y: &SimpleState, dy: &mut SimpleState) {
        let tm: AstroTime = self.start + x / 86400.0;
        *dy = SimpleState::zeros();
        for v in self.forces.iter() {
            *dy = *dy + v.ydot(&tm, y);
        }
        dy[0] += y[3];
        dy[1] += y[4];
        dy[2] += y[5];
    }
}

pub fn propagate(
    state: &SimpleState,
    start: &AstroTime,
    stop: &AstroTime,
    settings: &PropSettings,
) -> AstroResult<PropagationResult<SimpleState>> {
    // Gravity force recorder
    let g = &ForceGravity::init(start, stop, settings);
    let m = &ForceMoon::init(start, stop, settings);
    let s: &ForceSun = &ForceSun::init(start, stop, settings);
    // Propagation structure
    let p: Propagation = Propagation {
        start: *start,
        forces: vec![g, m, s],
    };

    // Duration to end of integration, in seconds
    let x_end: f64 = (*stop - *start) * 86400.0;

    // ODE stepper object
    let mut stepper: Dopri5<SimpleState, Propagation> = Dopri5::new(
        p,
        0.0,
        x_end,
        settings.dt_secs,
        *state,
        settings.rel_error,
        settings.abs_error,
    );

    // Perform the integration
    match stepper.integrate() {
        Ok(stats) => Ok(PropagationResult {
            time: stepper
                .x_out()
                .into_iter()
                .map(|x| *start + (x / 86400.0))
                .collect(),
            state: stepper.y_out().clone(),
            accepted_steps: stats.accepted_steps,
            rejected_steps: stats.rejected_steps,
            num_eval: stats.num_eval,
        }),
        Err(err) => {
            let s = err.to_string();
            return astroerr!("Propagation Error: {}", s);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::univ;

    #[test]
    fn test_propagate() {
        let starttime = AstroTime::from_datetime(2012, 3, 3, 0, 0, 0.0);
        let stoptime = starttime + 1.0;

        let mut state: SimpleState = SimpleState::zeros();
        state[0] = univ::GEO_R;
        state[4] = f64::sqrt(univ::MU_EARTH / univ::GEO_R);
        let mut settings = PropSettings::new();
        settings.dt_secs = 86400.0;
        let state2 = propagate(&state, &starttime, &stoptime, &settings);
        println!("state2 = {:?}", state2.unwrap().time);
    }
}
