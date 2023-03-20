use super::settings::PropSettings;

use crate::lpephem::{moon, sun};

use crate::astrotime::AstroTime;
use crate::utils::{astroerr, AstroResult};
use ode_solvers::dop853::Dop853;
//use ode_solvers::dopri5::Dopri5;
use std::vec::Vec;

use super::satproperties::SatProperties;
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

// T = 1 is simple state
// T = 42 is simple state + 6x6 state transition matrix
pub type StateType<const T: usize> = na::SVector<f64, T>;

pub type SimpleState = StateType<6>;

// Covariance State in includes
// 6x1 simple state and
// 6x6 state transition matrix
// for total of 42 elements
pub type CovState = StateType<42>;

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
    satprops: Option<&'a dyn SatProperties>,
}

impl<'a, const T: usize> ode_solvers::System<StateType<T>> for Propagation<'a> {
    fn system(&self, x: f64, y: &StateType<T>, dy: &mut StateType<T>) {
        let tm: AstroTime = self.start + x / 86400.0;

        // get GCRS position & velocity;
        let pos_gcrs: na::Vector3<f64> = y.fixed_slice::<3, 1>(0, 0).into();
        let vel_gcrs: na::Vector3<f64> = y.fixed_slice::<3, 1>(3, 0).into();

        // Get force of moon from interpolation table
        let moon_idx: f64 = (x / self.settings.moon_interp_dt_secs).floor();
        let moon_gcrs = linterp_idx(&self.moon_pos_gcrs_table, moon_idx).unwrap();

        // Get sun location & force of sun from interpolation table
        let sun_idx: f64 = x / self.settings.sun_interp_dt_secs;
        let sun_gcrs = linterp_idx(&self.sun_pos_gcrs_table, sun_idx).unwrap();

        // Get rotation from gcrf to itrf frame from interpolation table
        let grav_idx: usize = (x / self.settings.gravity_interp_dt_secs).floor() as usize;
        // t should be between 0 & 1
        let t = (x / self.settings.gravity_interp_dt_secs) - grav_idx as f64;
        let q1 = &self.qgcrs2itrf_table[grav_idx];
        let q2 = &self.qgcrs2itrf_table[grav_idx + 1];
        // Quaternion to go from inertial to terrestrial frame
        let qgcrs2itrf = q1.slerp(q2, t);

        // Position in ITRF coordinates
        let pos_itrf = qgcrs2itrf * pos_gcrs;

        // Propagating a "simple" 6-dof (position, velocity) state
        if T == 6 {
            // Gravity in the ITRF frame
            let gravity_itrf =
                crate::gravity::GRAVITY_JGM3.accel(&pos_itrf, self.settings.gravity_order as usize);
            // Gravity in the GCRS frame
            let accel_gravity = qgcrs2itrf.conjugate() * gravity_itrf;

            // Acceleration due to moon
            let accel_moon = point_gravity(&pos_gcrs, &moon_gcrs, crate::univ::MU_MOON);

            // Acceleration due to sun
            let accel_sun = point_gravity(&pos_gcrs, &sun_gcrs, crate::univ::MU_SUN);

            // Total acceleration (neglecting solar pressure & drag)
            let mut accel = accel_gravity + accel_sun + accel_moon;

            // Add solar pressure & drag if that is defined in satellite properties
            if let Some(props) = self.satprops {
                let ss = y.fixed_slice::<6, 1>(0, 0);

                // Compute solar pressure
                let solarpressure = -props.cr_a_over_m(&tm, &ss.into()) * 4.56e-6 * sun_gcrs
                    / sun_gcrs.norm().powf(3.0);
                accel += solarpressure;

                // Compute drag
                let cd_a_over_m = props.cd_a_over_m(&tm, &ss.into());
                if cd_a_over_m > 1e-6 {
                    let itrf = crate::ITRFCoord::from(pos_itrf.as_slice());
                    let (density, _temperature) = crate::nrlmsise::nrlmsise(
                        itrf.hae() / 1.0e3,
                        Some(itrf.latitude_rad()),
                        Some(itrf.longitude_rad()),
                        Some(tm),
                        self.settings.use_spaceweather,
                    );
                    const OMEGA_EARTH: na::Vector3<f64> =
                        na::vector![0.0, 0.0, crate::univ::OMEGA_EARTH];

                    // See Vallado section 3.7.2: Velocity & Acceleration Transformations
                    let vel_itrf = qgcrs2itrf * vel_gcrs - OMEGA_EARTH.cross(&pos_itrf);
                    let dragpressure_itrf =
                        -0.5 * cd_a_over_m * density * vel_itrf * vel_itrf.norm();
                    let dragpressure_gcrs = qgcrs2itrf.conjugate()
                        * (dragpressure_itrf
                            + OMEGA_EARTH.cross(&OMEGA_EARTH.cross(&pos_itrf))
                            + 2.0 * OMEGA_EARTH.cross(&vel_itrf));
                    accel += dragpressure_gcrs;
                }
            } // end of handling drag & solarpressure

            // change in position is velocity
            dy.fixed_slice_mut::<3, 1>(0, 0).copy_from(&vel_gcrs);

            // Change in velocity is acceleration
            dy.fixed_slice_mut::<3, 1>(3, 0).copy_from(&accel);
        }
    }
}

pub fn propagate<const T: usize>(
    state: &StateType<T>,
    start: &AstroTime,
    stop: &AstroTime,
    settings: &PropSettings,
    satprops: Option<&dyn SatProperties>,
) -> AstroResult<PropagationResult<StateType<T>>> {
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
        satprops: satprops,
    };

    // Duration to end of integration, in seconds
    let x_end: f64 = (*stop - *start) * 86400.0;

    // ODE stepper object
    let mut stepper: Dop853<StateType<T>, Propagation> = Dop853::<StateType<T>, Propagation>::new(
        p,
        0.0,
        x_end,
        settings.dt_secs,
        state.clone(),
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
    fn test_nalgebra() {
        let v = na::SVector::<f64, 42>::from_iterator((0..42).map(|x| x as f64));
        let v2 = v.reshape_generic(na::Const::<6>, na::Const::<7>);
        println!("v = {}", v.transpose());
        println!("v2 = {v2}");
    }

    #[test]
    fn test_propagate() {
        let starttime = AstroTime::from_datetime(2015, 3, 20, 0, 0, 0.0);
        let stoptime = starttime + 1.0; // - (1.0 / 365.25);

        let mut state: SimpleState = SimpleState::zeros();
        state[0] = univ::GEO_R;
        state[4] = (univ::MU_EARTH / univ::GEO_R).sqrt();
        let mut settings = PropSettings::new();
        settings.dt_secs = 1.0;
        settings.gravity_order = 6;
        let state2 = propagate(&state, &starttime, &stoptime, &settings, None).unwrap();

        println!("start time = {}", starttime);
        println!("start = {}", state);
        println!("end time = {}", &state2.time.last().unwrap());
        println!("end = {}", state2.state.last().unwrap());
        let pg1 = state2.state.first().unwrap();
        let pg2 = state2.state.last().unwrap();

        let q1 = crate::frametransform::qgcrf2itrf(&state2.time.first().unwrap());
        let q2 = crate::frametransform::qgcrf2itrf(&state2.time.last().unwrap());
        let p1itrf = q1 * na::vector![pg1[0], pg1[1], pg1[2]];
        let p2itrf = q2 * na::vector![pg2[0], pg2[1], pg2[2]];

        let p1 = crate::itrfcoord::ITRFCoord::from_vec([p1itrf[0], p1itrf[1], p1itrf[2]]);
        let p2 = crate::itrfcoord::ITRFCoord::from_vec([p2itrf[0], p2itrf[1], p2itrf[2]]);

        println!("p1 = {}", p1);
        println!("p2 = {}", p2);
        println!("pdiff = {}", p1itrf - p2itrf);
        println!("rejected = {}", state2.rejected_steps);
        println!("accepted = {}", state2.accepted_steps);
    }
}
