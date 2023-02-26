use super::settings::PropSettings;
use super::types::*;

use crate::lpephem::{moon, sun};

use crate::astrotime::AstroTime;
use crate::utils::{astroerr, AstroResult};
use ode_solvers::dop853::Dop853;
use ode_solvers::dopri5::Dopri5;
use std::vec::Vec;

use super::utils::linterp_idx;

use nalgebra as na;

#[derive(Debug)]
pub struct PropagationResult<T> {
    pub time: Vec<AstroTime>,
    pub state: Vec<T>,
    pub accepted_steps: u32,
    pub rejected_steps: u32,
    pub num_eval: u32,
}

pub type CovState = na::Matrix<f64, na::Const<6>, na::Const<7>, na::ArrayStorage<f64, 6, 7>>;

// T = 1 is simple state
// T = 7 is simple state + 6x6 state transition matrix
pub type StateType<const T: usize> =
    na::Matrix<f64, na::Const<6>, na::Const<T>, na::ArrayStorage<f64, 6, T>>;

// Equation 3.37 in Montenbruck & Gill
fn point_gravity(
    r: &na::Vector3<f64>, // object
    s: &na::Vector3<f64>, // distant attractor
    mu: f64,
) -> na::Vector3<f64> {
    let rs = r - s;
    let rsnorm2 = rs.norm_squared();
    let rsnorm = rsnorm2.sqrt();
    let snorm2 = s.norm_squared();
    let snorm = snorm2.sqrt();
    -mu * ((rs / (rsnorm * rsnorm2)) + (s / (snorm * snorm2)))
}

struct Propagation<'a> {
    start: AstroTime,
    settings: &'a PropSettings,
    sun_pos_gcrs_table: Vec<na::Vector3<f64>>,
    moon_pos_gcrs_table: Vec<na::Vector3<f64>>,
    qgcrs2itrf_table: Vec<na::UnitQuaternion<f64>>,
}

impl<'a, const T: usize> ode_solvers::System<StateType<T>> for Propagation<'a> {
    fn system(&self, x: f64, y: &StateType<T>, dy: &mut StateType<T>) {
        //let tm: AstroTime = self.start + x / 86400.0;
        //*dy = StateType::<T>::zeros();

        // get GCRS position & velocity;
        let pos_gcrs: na::Vector3<f64> = y.fixed_slice::<3, 1>(0, 0).into();
        let vel_gcrs: na::Vector3<f64> = y.fixed_slice::<3, 1>(3, 0).into();

        // Get force of moon from interpolation table
        let moon_idx: f64 = (x / self.settings.moon_interp_dt_secs).floor();

        let moon_gcrs = linterp_idx(&self.moon_pos_gcrs_table, moon_idx).unwrap();
        let accel_moon = point_gravity(&pos_gcrs, &moon_gcrs, crate::univ::MU_MOON);

        // Get sun location & force of sun from interpolation table
        let sun_idx: f64 = x / self.settings.sun_interp_dt_secs;
        let sun_gcrs = linterp_idx(&self.sun_pos_gcrs_table, sun_idx).unwrap();
        let accel_sun = point_gravity(&pos_gcrs, &sun_gcrs, crate::univ::MU_SUN);

        // Get rotation from gcrf to itrf frame from interpolation table
        let grav_idx: usize = (x / self.settings.gravity_interp_dt_secs).floor() as usize;
        let t = (x / self.settings.gravity_interp_dt_secs) - grav_idx as f64;
        let q1 = &self.qgcrs2itrf_table[grav_idx];
        let q2 = &self.qgcrs2itrf_table[grav_idx + 1];
        let qgcrs2itrf = q1.slerp(q2, t);

        let pos_itrf = qgcrs2itrf * pos_gcrs;
        let gravity_itrf =
            crate::gravity::GRAVITY_JGM3.accel(&pos_itrf, self.settings.gravity_order as usize);
        let accel_gravity = qgcrs2itrf.conjugate() * gravity_itrf;

        // change in position is velocity
        dy.fixed_slice_mut::<3, 1>(0, 0).copy_from(&vel_gcrs);

        // Change in velocity is acceleration
        let accel = accel_gravity + accel_sun + accel_moon;
        //let accel = -crate::univ::MU_EARTH * pos_gcrs / pos_gcrs.norm().powf(3.0);

        dy.fixed_slice_mut::<3, 1>(3, 0).copy_from(&accel);
    }
}

pub fn propagate(
    state: &SimpleState,
    start: &AstroTime,
    stop: &AstroTime,
    settings: &PropSettings,
) -> AstroResult<PropagationResult<SimpleState>> {
    // Propagation structure
    let duration_days: f64 = *stop - *start;
    let duration_secs: f64 = duration_days * 86400.0;

    let p: Propagation = Propagation {
        start: *start,
        settings: settings,

        // Sun positions for interpolation
        sun_pos_gcrs_table: {
            let ntimes: u32 = 2 + (duration_secs / settings.sun_interp_dt_secs).ceil() as u32;
            (0..ntimes)
                .into_iter()
                .map(|x| {
                    let tm: AstroTime = *start + x as f64 / 86400.0;
                    sun::pos_gcrf(&tm)
                })
                .collect()
        },

        // Moon positions for interpolation
        moon_pos_gcrs_table: {
            let ntimes: u32 = 2 + (duration_secs / settings.moon_interp_dt_secs).ceil() as u32;
            (0..ntimes)
                .into_iter()
                .map(|x| {
                    let tm: AstroTime = *start + x as f64 / 86400.0;
                    moon::pos_gcrf(&tm)
                })
                .collect()
        },

        // qGCRS2ITRf
        qgcrs2itrf_table: {
            let ntimes: u32 = 2 + (duration_secs / settings.gravity_interp_dt_secs).ceil() as u32;
            (0..ntimes)
                .into_iter()
                .map(|x| {
                    let tm: AstroTime = *start + x as f64 / 86400.0;
                    crate::frametransform::qgcrf2itrf_approx(&tm)
                })
                .collect()
        },
    };

    // Duration to end of integration, in seconds
    let x_end: f64 = (*stop - *start) * 86400.0;
    println!("x_end = {x_end}");

    // ODE stepper object
    let mut stepper: Dop853<StateType<1>, Propagation> = Dop853::new(
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
        state[4] = (univ::MU_EARTH / univ::GEO_R).sqrt();
        let mut settings = PropSettings::new();
        settings.dt_secs = 1.0;
        settings.gravity_order = 6;
        let state2 = propagate(&state, &starttime, &stoptime, &settings).unwrap();

        println!("state star ttime = {:?}", starttime);
        println!("state start = {:?}", state);
        println!("state2 end time = {:?}", &state2.time.last().unwrap());
        println!("dt = {}", *state2.time.last().unwrap() - starttime);
        println!("state2 end state = {:?}", state2.state.last().unwrap());
        let pg1 = state2.state.first().unwrap();
        let pg2 = state2.state.last().unwrap();

        let q1 = crate::frametransform::qgcrf2itrf(&state2.time.first().unwrap());
        let q2 = crate::frametransform::qgcrf2itrf(&state2.time.last().unwrap());
        let p1itrf = q1 * na::vector![pg1[0], pg1[1], pg1[2]];
        let p2itrf = q2 * na::vector![pg2[0], pg2[1], pg2[2]];

        let p1 = crate::itrfcoord::ITRFCoord::from_vec([p1itrf[0], p1itrf[1], p1itrf[2]]);
        let p2 = crate::itrfcoord::ITRFCoord::from_vec([p2itrf[0], p2itrf[1], p2itrf[2]]);
        println!("q1 = {q1}");
        println!("q2 = {q2}");
        println!("p1 = {}", p1);
        println!("p2 = {}", p2);
        println!("pdiff = {:?}", p1itrf - p2itrf);
    }
}
