use super::settings::PropSettings;

use crate::astrotime::AstroTime;
use crate::ode;
use crate::ode::RKAdaptive;

use crate::utils::AstroResult;

use num_traits::identities::Zero;

use std::vec::Vec;

use super::utils::linterp_idx;
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

#[derive(Debug, Error)]
pub enum PropagationError {
    #[error("Invalid number of columns: {c}")]
    InvalidStateColumns { c: usize },
}

struct Propagation<'a, const C: usize> {
    start: AstroTime,
    settings: &'a PropSettings,
    sun_pos_gcrs_table: Vec<na::Vector3<f64>>,
    moon_pos_gcrs_table: Vec<na::Vector3<f64>>,
    qgcrs2itrf_table: Vec<na::UnitQuaternion<f64>>,
    satprops: Option<&'a dyn SatProperties>,
}

impl<'a, const C: usize> ode::ODESystem for Propagation<'a, C> {
    type Output = StateType<C>;

    fn ydot(&mut self, x: f64, y: &Self::Output) -> ode::ODEResult<Self::Output> {
        let tm: AstroTime = self.start + x / 86400.0;

        // get GCRS position & velocity;
        let pos_gcrs: na::Vector3<f64> = y.fixed_view::<3, 1>(0, 0).into();
        let vel_gcrs: na::Vector3<f64> = y.fixed_view::<3, 1>(3, 0).into();

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
        if C == 1 {
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
            let mut accel = accel_gravity - accel_sun - accel_moon;

            // Add solar pressure & drag if that is defined in satellite properties
            if let Some(props) = self.satprops {
                let ss = y.fixed_view::<6, 1>(0, 0);

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

            let mut dy = Self::Output::zero();
            // change in position is velocity
            dy.fixed_view_mut::<3, 1>(0, 0).copy_from(&vel_gcrs);

            // Change in velocity is acceleration
            dy.fixed_view_mut::<3, 1>(3, 0).copy_from(&accel);

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
    settings: &PropSettings,
    satprops: Option<&dyn SatProperties>,
) -> AstroResult<crate::ode::ODESolution<StateType<C>>> {
    // Propagation structure
    let duration_days: f64 = *stop - *start;
    let duration_secs: f64 = duration_days * 86400.0;

    let mut p: Propagation<C> = Propagation::<C> {
        start: *start,
        settings: settings,

        // Sun positions for interpolation
        sun_pos_gcrs_table: {
            let ntimes: u32 = 2 + (duration_secs / settings.sun_interp_dt_secs).ceil() as u32;
            (0..ntimes)
                .into_iter()
                .map(|x| {
                    let tm: AstroTime = *start + x as f64 * settings.sun_interp_dt_secs / 86400.0;
                    crate::jplephem::geocentric_pos(crate::SolarSystem::SUN, &tm).unwrap()
                })
                .collect()
        },

        // Moon positions for interpolation
        moon_pos_gcrs_table: {
            let ntimes: u32 = 2 + (duration_secs / settings.moon_interp_dt_secs).ceil() as u32;
            (0..ntimes)
                .into_iter()
                .map(|x| {
                    let tm: AstroTime = *start + x as f64 * settings.moon_interp_dt_secs / 86400.0;
                    crate::jplephem::geocentric_pos(crate::SolarSystem::MOON, &tm).unwrap()
                })
                .collect()
        },

        // qGCRS2ITRf
        qgcrs2itrf_table: {
            let ntimes: u32 = 2 + (duration_secs / settings.gravity_interp_dt_secs).ceil() as u32;
            (0..ntimes)
                .into_iter()
                .map(|x| {
                    let tm: AstroTime =
                        *start + x as f64 * settings.gravity_interp_dt_secs / 86400.0;
                    crate::frametransform::qgcrf2itrf_approx(&tm)
                })
                .collect()
        },
        satprops: satprops,
    };

    // Duration to end of integration, in seconds
    let x_end: f64 = (*stop - *start) * 86400.0;

    let mut settings = crate::ode::RKAdaptiveSettings::default();
    settings.abserror = 1.0e-2;
    settings.relerror = 1.0e-2;
    println!("integrating");
    let res = crate::ode::RKV65::integrate(0.0, x_end, state, &mut p, &settings)?;

    Ok(res)
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

        let mut settings = PropSettings::new();
        settings.dt_secs = 1.0;
        settings.gravity_order = 6;
        settings.gravity_interp_dt_secs = 60.0;
        settings.moon_interp_dt_secs = 60.0;
        settings.sun_interp_dt_secs = 60.0;

        println!("running");
        let res = propagate(&state, &starttime, &stoptime, &settings, None)?;
        println!("res = {:?}", res);
        Ok(())
    }
}
