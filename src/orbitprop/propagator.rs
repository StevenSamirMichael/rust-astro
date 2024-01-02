use super::settings::PropSettings;

use crate::astrotime::AstroTime;
use crate::earthgravity;
use crate::frametransform;
use crate::jplephem;
use crate::lpephem;
use crate::ode;
use crate::ode::RKAdaptive;
use crate::SolarSystem;

use crate::utils::AstroResult;

use num_traits::identities::Zero;

use std::vec::Vec;

use crate::orbitprop::SatProperties;
use crate::univ;

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

// Equation 3.37 in Montenbruck & Gill
fn point_gravity(
    r: &na::Vector3<f64>, // object
    s: &na::Vector3<f64>, // distant attractor
    mu: f64,
) -> na::Vector3<f64> {
    let sr = s - r;
    let srnorm2 = sr.norm_squared();
    let srnorm = srnorm2.sqrt();
    let snorm2 = s.norm_squared();
    let snorm = snorm2.sqrt();
    mu * (sr / (srnorm * srnorm2) - s / (snorm * snorm2))
}

// Return tuple with point gravity force and
// point gravity partial (da/dr)
// Equation 3.37 in Montenbruck & Gill for point gravity
// Equation 7.75 in Montenbruck & Gill for partials

fn point_gravity_and_partials(
    r: &na::Vector3<f64>, // object
    s: &na::Vector3<f64>, // distant attractur
    mu: f64,
) -> (na::Vector3<f64>, na::Matrix3<f64>) {
    let rs = r - s;
    let rsnorm2 = rs.norm_squared();
    let rsnorm = rsnorm2.sqrt();
    let snorm2 = s.norm_squared();
    let snorm = snorm2.sqrt();
    (
        -mu * (rs / (rsnorm * rsnorm2) + s / (snorm * snorm2)),
        -mu * (na::Matrix3::<f64>::identity() / (rsnorm2 * rsnorm)
            - 3.0 * rs * rs.transpose() / (rsnorm2 * rsnorm2 * rsnorm)),
    )
}

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

impl<'a, const C: usize> ode::ODESystem for Propagation<'a, C> {
    type Output = StateType<C>;

    fn ydot(&mut self, x: f64, y: &Self::Output) -> ode::ODEResult<Self::Output> {
        let tm: AstroTime = self.start + x / 86400.0;
        //println!("tm = {} :: {} :: x = {}", self.start, tm, x);

        // get GCRS position & velocity;
        let pos_gcrf: na::Vector3<f64> = y.fixed_view::<3, 1>(0, 0).into();
        let vel_gcrf: na::Vector3<f64> = y.fixed_view::<3, 1>(3, 0).into();

        // Moon position
        let (sun_gcrf, moon_gcrf) = match self.settings.use_jplephem {
            true => (
                jplephem::geocentric_pos(SolarSystem::SUN, &tm)?,
                jplephem::geocentric_pos(SolarSystem::MOON, &tm)?,
            ),
            false => (lpephem::sun::pos_gcrf(&tm), lpephem::moon::pos_gcrf(&tm)),
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
            accel += point_gravity(&pos_gcrf, &moon_gcrf, crate::univ::MU_MOON);

            // Acceleration due to sun
            accel += point_gravity(&pos_gcrf, &sun_gcrf, crate::univ::MU_SUN);

            // Add solar pressure & drag if that is defined in satellite properties
            if let Some(props) = self.satprops {
                let ss = y.fixed_view::<6, 1>(0, 0);

                // Compute solar pressure
                let solarpressure =
                    -props.cr_a_over_m(&tm, &ss.into()) * 4.56e-6 * sun_gcrf / sun_gcrf.norm();
                accel += solarpressure;

                // Compute drag
                if pos_gcrf.norm() < 700.0e3 + crate::univ::EARTH_RADIUS {
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
                        let vel_itrf = qgcrf2itrf * vel_gcrf - OMEGA_EARTH.cross(&pos_itrf);

                        let dragpressure_itrf =
                            -0.5 * cd_a_over_m * density * vel_itrf * vel_itrf.norm();
                        let dragpressure_gcrs = qgcrf2itrf.conjugate()
                            * (dragpressure_itrf
                                + OMEGA_EARTH.cross(&OMEGA_EARTH.cross(&pos_itrf))
                                + 2.0 * OMEGA_EARTH.cross(&vel_itrf));
                        accel += dragpressure_gcrs;
                    }
                }
            } // end of handling drag & solarpressure

            let mut dy = Self::Output::zero();
            // change in position is velocity
            dy.fixed_view_mut::<3, 1>(0, 0).copy_from(&vel_gcrf);

            // Change in velocity is acceleration
            dy.fixed_view_mut::<3, 1>(3, 0).copy_from(&accel);

            Ok(dy)
        } else if C == 7 {
            let (gravity_accel, gravity_partials) = earthgravity::jgm3()
                .accel_and_partials(&pos_itrf, self.settings.gravity_order as usize);
            let (sun_accel, sun_partials) =
                point_gravity_and_partials(&pos_gcrf, &sun_gcrf, univ::MU_SUN);
            let (moon_accel, moon_partials) =
                point_gravity_and_partials(&pos_gcrf, &moon_gcrf, univ::MU_MOON);

            let accel = qitrf2gcrf * gravity_accel + sun_accel + moon_accel;

            // Equation 7.42 in Montenbruck & Gill
            let mut dfdy: StateType<6> = StateType::<6>::zeros();
            dfdy.fixed_view_mut::<3, 3>(0, 3)
                .copy_from(&na::Matrix3::<f64>::identity());

            let ritrf2gcrf = qitrf2gcrf.to_rotation_matrix();
            let dadr = ritrf2gcrf * gravity_partials * ritrf2gcrf.transpose()
                + sun_partials
                + moon_partials;

            dfdy.fixed_view_mut::<3, 3>(3, 0).copy_from(&dadr);

            let dphi = dfdy * y.fixed_view::<6, 6>(0, 1);

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
    let duration_days: f64 = (*stop - *start).abs();
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
    let x_end: f64 = (*stop - *start) * 86400.0;

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
    use crate::univ;

    #[test]
    fn test_propagate() -> AstroResult<()> {
        let starttime = AstroTime::from_datetime(2015, 3, 20, 0, 0, 0.0);
        let stoptime = starttime + 1.0;

        let mut state: SimpleState = SimpleState::zeros();
        state[0] = univ::GEO_R;
        state[4] = (univ::MU_EARTH / univ::GEO_R).sqrt();

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
        state[0] = univ::GEO_R;
        state[4] = (univ::MU_EARTH / univ::GEO_R).sqrt();
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
