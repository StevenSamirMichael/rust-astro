use super::drag::{drag_and_partials, drag_force};
use super::point_gravity::{point_gravity, point_gravity_and_partials};
use super::settings::PropSettings;

use crate::astrotime::AstroTime;
use crate::earthgravity;
use crate::frametransform;
use crate::jplephem;
use crate::lpephem;
use crate::ode;
use crate::ode::RKAdaptive;
use crate::Duration;
use crate::SolarSystem;
use lpephem::sun::shadowfunc;

use crate::utils::AstroResult;

use num_traits::identities::Zero;

use std::vec::Vec;

use crate::consts;
use crate::orbitprop::SatProperties;

use thiserror::Error;

use nalgebra as na;

#[derive(Debug)]
pub struct PropagationResult<T> {
    pub time: Vec<AstroTime>,
    pub state: Vec<T>,
    pub accepted_steps: u32,
    pub rejected_steps: u32,
    pub num_eval: u32,
}

pub type StateType<const C: usize> = na::SMatrix<f64, 6, C>;
pub type SimpleState = StateType<1>;

// Covariance State in includes
pub type CovState = StateType<7>;

#[derive(Debug, Error)]
pub enum PropagationError {
    #[error("Invalid number of columns: {c}")]
    InvalidStateColumns { c: usize },
}

struct Propagation<'a, const C: usize> {
    start: AstroTime,
    settings: &'a PropSettings,
    qgcrs2itrf_table: Vec<na::UnitQuaternion<f64>>,
    satprops: Option<&'a dyn SatProperties>,
}

//
// This actually implements the force model that is used to
// integrate the ODE to get position and velocity
//
// State is position and velocity
// Force is computed and integrated to get velocity
// Velocity is integrated to get position
//
// If C=7, a 6x6 state transition matrix is appended as additional
// colums, making the integrated "state" a 6x7 matrix
// The state transition matrix can be used to propagate covariances
//
// See Montenbruk & Gill for details (Chapter 7)
//
impl<'a, const C: usize> ode::ODESystem for Propagation<'a, C> {
    type Output = StateType<C>;

    fn ydot(&mut self, x: f64, y: &Self::Output) -> ode::ODEResult<Self::Output> {
        // The time variable in the ODE is in seconds
        let time: AstroTime = self.start + Duration::Seconds(x);

        // get GCRS position & velocity;
        let pos_gcrf: na::Vector3<f64> = y.fixed_view::<3, 1>(0, 0).into();
        let vel_gcrf: na::Vector3<f64> = y.fixed_view::<3, 1>(3, 0).into();

        // Moon position
        let (sun_gcrf, moon_gcrf) = match self.settings.use_jplephem {
            true => (
                jplephem::geocentric_pos(SolarSystem::SUN, &time)?,
                jplephem::geocentric_pos(SolarSystem::MOON, &time)?,
            ),
            false => (
                lpephem::sun::pos_gcrf(&time),
                lpephem::moon::pos_gcrf(&time),
            ),
        };

        // Get rotation from gcrf to itrf frame from interpolation table
        let grav_idx: f64 = x.abs() / self.settings.gravity_interp_dt_secs;
        // t should be between 0 & 1
        let t = grav_idx - grav_idx.floor();

        let q1 = &self.qgcrs2itrf_table[grav_idx as usize];
        let q2 = &self.qgcrs2itrf_table[grav_idx as usize + 1];
        // Quaternion to go from inertial to terrestrial frame
        let qgcrf2itrf = q1.slerp(q2, t);
        //let qgcrf2itrf = frametransform::qgcrf2itrf_approx(&tm);
        let qitrf2gcrf = qgcrf2itrf.conjugate();

        // Position in ITRF coordinates
        let pos_itrf = qgcrf2itrf * pos_gcrf;

        // Propagating a "simple" 6-dof (position, velocity) state
        if C == 1 {
            let mut accel = na::Vector3::<f64>::zeros();

            // Gravity in the ITRF frame
            let gravity_itrf =
                earthgravity::jgm3().accel(&pos_itrf, self.settings.gravity_order as usize);

            // Gravity in the GCRS frame
            accel += qitrf2gcrf * gravity_itrf;

            // Acceleration due to moon
            accel += point_gravity(&pos_gcrf, &moon_gcrf, crate::consts::MU_MOON);

            // Acceleration due to sun
            accel += point_gravity(&pos_gcrf, &sun_gcrf, crate::consts::MU_SUN);

            // Add solar pressure & drag if that is defined in satellite properties
            if let Some(props) = self.satprops {
                let ss = y.fixed_view::<6, 1>(0, 0);

                // Compute solar pressure
                let solarpressure = -shadowfunc(&sun_gcrf, &pos_gcrf)
                    * props.cr_a_over_m(&time, &ss.into())
                    * 4.56e-6
                    * sun_gcrf
                    / sun_gcrf.norm();
                accel += solarpressure;

                // Compute drag
                if pos_gcrf.norm() < 700.0e3 + crate::consts::EARTH_RADIUS {
                    let cd_a_over_m = props.cd_a_over_m(&time, &ss.into());

                    if cd_a_over_m > 1e-6 {
                        accel += drag_force(
                            &pos_gcrf,
                            &pos_itrf,
                            &vel_gcrf,
                            &time,
                            cd_a_over_m,
                            self.settings.use_spaceweather,
                        );
                    }
                }
            } // end of handling drag & solarpressure

            let mut dy = Self::Output::zero();
            // change in position is velocity
            dy.fixed_view_mut::<3, 1>(0, 0).copy_from(&vel_gcrf);

            // Change in velocity is acceleration
            dy.fixed_view_mut::<3, 1>(3, 0).copy_from(&accel);

            Ok(dy)
        }
        // If C==7, we are also integrating the state transition matrix
        else if C == 7 {
            // For state transition matrix, we need to compute force partials with respect to position
            // (for all forces but drag, partial with respect to velocity are zero)
            let (gravity_accel, gravity_partials) = earthgravity::jgm3()
                .accel_and_partials(&pos_itrf, self.settings.gravity_order as usize);
            let (sun_accel, sun_partials) =
                point_gravity_and_partials(&pos_gcrf, &sun_gcrf, consts::MU_SUN);
            let (moon_accel, moon_partials) =
                point_gravity_and_partials(&pos_gcrf, &moon_gcrf, consts::MU_MOON);

            let mut accel = qitrf2gcrf * gravity_accel + sun_accel + moon_accel;

            // Equation 7.42 in Montenbruck & Gill
            let mut dfdy: StateType<6> = StateType::<6>::zeros();
            dfdy.fixed_view_mut::<3, 3>(0, 3)
                .copy_from(&na::Matrix3::<f64>::identity());

            let ritrf2gcrf = qitrf2gcrf.to_rotation_matrix();
            // Sum partials with respect to position for gravity, sun, and moon
            // Note: gravity partials need to be rotated into the gcrf frame from itrf
            let mut dadr = ritrf2gcrf * gravity_partials * ritrf2gcrf.transpose()
                + sun_partials
                + moon_partials;

            // Handle satellite properties for drag and radiation pressure
            if let Some(props) = self.satprops {
                // Satellite state as 6-element position, velcoity matrix
                // used to query cd_a_over_m
                let ss = y.fixed_view::<6, 1>(0, 0);

                // Compute solar pressure
                // Partials for this are very small since the sun is very very far away, changes in
                // satellite position don't change radiaion pressure much, so we will ignore...
                let solarpressure = -shadowfunc(&sun_gcrf, &pos_gcrf)
                    * props.cr_a_over_m(&time, &ss.into())
                    * 4.56e-6
                    * sun_gcrf
                    / sun_gcrf.norm();
                accel += solarpressure;

                // We know drag is negligible above 700 km, so ignore if this is the case
                if pos_gcrf.norm() < 700.0e3 + crate::consts::EARTH_RADIUS {
                    let cd_a_over_m = props.cd_a_over_m(&time, &ss.into());
                    if cd_a_over_m > 1e-6 {
                        let (drag_accel, ddragaccel_dr, ddragaccel_dv) = drag_and_partials(
                            &pos_gcrf,
                            &qgcrf2itrf,
                            &vel_gcrf,
                            &time,
                            cd_a_over_m,
                            self.settings.use_spaceweather,
                        );

                        // Add acceleration from drag to accel vector
                        accel += drag_accel;

                        // Add drag partials with respect to position to
                        // daccel dr
                        dadr += ddragaccel_dr;

                        // Drag is the only force term that produces a finite partial with respect
                        // to velocity, so copy it directly into dfdy here.
                        dfdy.fixed_view_mut::<3, 3>(3, 3).copy_from(&ddragaccel_dv);
                    }
                }
            }
            dfdy.fixed_view_mut::<3, 3>(3, 0).copy_from(&dadr);

            // Derivative of state transition matrix is dfdy * state transition matrix
            let dphi: na::Matrix<f64, na::Const<6>, na::Const<6>, na::ArrayStorage<f64, 6, 6>> =
                dfdy * y.fixed_view::<6, 6>(0, 1);

            let mut dy = Self::Output::zero();
            dy.fixed_view_mut::<3, 1>(0, 0).copy_from(&vel_gcrf);
            dy.fixed_view_mut::<3, 1>(3, 0).copy_from(&accel);
            dy.fixed_view_mut::<6, 6>(0, 1).copy_from(&dphi);
            Ok(dy)
        } else {
            Err(Box::new(PropagationError::InvalidStateColumns { c: C }))
        }
    }
}

pub fn propagate<const C: usize>(
    state: &StateType<C>,
    start: &AstroTime,
    stop: &AstroTime,
    step_seconds: Option<f64>,
    settings: &PropSettings,
    satprops: Option<&dyn SatProperties>,
) -> AstroResult<PropagationResult<StateType<C>>> {
    // Propagation structure
    let duration_days: f64 = (*stop - *start).days().abs();
    let duration_secs: f64 = duration_days * 86400.0;

    // Direction of time
    let tdir = match *stop > *start {
        true => 1.0,
        false => -1.0,
    };

    let mut p: Propagation<C> = Propagation::<C> {
        start: *start,
        settings: settings,

        // qGCRS2ITRf
        qgcrs2itrf_table: {
            let ntimes: u32 = 2 + (duration_secs / settings.gravity_interp_dt_secs).ceil() as u32;
            (0..ntimes)
                .into_iter()
                .map(|x| {
                    let tm: AstroTime =
                        *start + x as f64 * tdir * settings.gravity_interp_dt_secs / 86400.0; // Use high-precision transform
                    frametransform::qgcrf2itrf(&tm)
                })
                .collect()
        },
        satprops: satprops,
    };

    // Duration to end of integration, in seconds
    let x_end: f64 = (*stop - *start).days() * 86400.0;

    let mut odesettings = crate::ode::RKAdaptiveSettings::default();
    odesettings.abserror = settings.abs_error;
    odesettings.relerror = settings.rel_error;

    match step_seconds {
        None => {
            odesettings.dense_output = false;
            // If no interpolation, run a different integrator with same order but skipping
            // stages (& this computation time) that are only used for interpolation
            let res = crate::ode::solvers::RKV98NoInterp::integrate(
                0.0,
                x_end,
                state,
                &mut p,
                &odesettings,
            )?;
            Ok(PropagationResult {
                time: vec![stop.clone()],
                state: vec![res.y],
                accepted_steps: res.naccept as u32,
                rejected_steps: res.nreject as u32,
                num_eval: res.nevals as u32,
            })
        }
        Some(dx) => {
            odesettings.dense_output = true;
            let (res, interp) = crate::ode::solvers::RKV98::integrate_dense(
                0.0,
                x_end,
                dx,
                state,
                &mut p,
                &odesettings,
            )?;
            Ok(PropagationResult {
                time: interp.x.iter().map(|x| *start + *x / 86400.0).collect(),
                state: interp.y,
                accepted_steps: res.naccept as u32,
                rejected_steps: res.nreject as u32,
                num_eval: res.nevals as u32,
            })
        }
    }

    //Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts;

    #[test]
    fn test_propagate() -> AstroResult<()> {
        let starttime = AstroTime::from_datetime(2015, 3, 20, 0, 0, 0.0);
        let stoptime = starttime + 1.0;

        let mut state: SimpleState = SimpleState::zeros();
        state[0] = consts::GEO_R;
        state[4] = (consts::MU_EARTH / consts::GEO_R).sqrt();

        let mut settings = PropSettings::default();
        settings.abs_error = 1.0e-9;
        settings.rel_error = 1.0e-14;
        settings.gravity_order = 4;
        settings.gravity_interp_dt_secs = 300.0;
        settings.use_jplephem = false;

        println!("state0 = {}", state.transpose());
        println!("running");
        let res = propagate(&state, &starttime, &stoptime, None, &settings, None)?;
        println!("res = {:?}", res);
        println!("time = {}", res.time.last().unwrap());

        // Try to propagate back to original time
        let res2 = propagate(&res.state[0], &stoptime, &starttime, None, &settings, None);
        println!("res2 = {:?}", res2);

        Ok(())
    }

    #[test]
    fn test_state_transition() -> AstroResult<()> {
        let starttime = AstroTime::from_datetime(2015, 3, 20, 0, 0, 0.0);
        let stoptime = starttime + 1.0;

        let mut state: CovState = CovState::zeros();
        state[0] = consts::GEO_R;
        state[4] = (consts::MU_EARTH / consts::GEO_R).sqrt();
        state
            .fixed_view_mut::<6, 6>(0, 1)
            .copy_from(&na::Matrix6::<f64>::identity());
        println!("state 0 = {}", state);

        let mut settings = PropSettings::default();
        settings.abs_error = 1.0e-9;
        settings.rel_error = 1.0e-14;
        settings.gravity_order = 4;
        settings.gravity_interp_dt_secs = 300.0;
        settings.use_jplephem = false;

        println!("running");
        let res = propagate(&state, &starttime, &stoptime, None, &settings, None)?;
        println!("res = {:?}", res);
        println!("time = {}", res.time.last().unwrap());
        println!("state end = {}", res.state[0]);

        let res2 = propagate(&res.state[0], &stoptime, &starttime, None, &settings, None)?;
        println!("res2 state = {}", res2.state[0]);

        Ok(())
    }
}
