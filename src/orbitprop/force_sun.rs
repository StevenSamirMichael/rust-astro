use super::forceterm::*;
use super::point_gravity::point_gravity;
use super::types::*;
use super::utils::linterp_idx;
use super::PropSettings;

use nalgebra as na;

use crate::lpephem::sun;
use crate::univ;
use crate::AstroTime;

pub struct ForceSun {
    start_time: AstroTime,
    sun_interp_dt_secs: f64,
    pos_gcrf: Vec<na::Vector3<f64>>,
}

const SECONDS_PER_DAY: f64 = 86400.0;

impl ForceTerm<SimpleState> for ForceSun {
    fn ydot(&self, time: &AstroTime, state: &SimpleState) -> SimpleState {
        let float_idx: f64 = (*time - self.start_time) / self.sun_interp_dt_secs;
        let sun_pos = linterp_idx(&self.pos_gcrf, float_idx).unwrap();
        let pos: Vec3 = Vec3::from_row_slice(state.index((0..3, ..)).as_slice());
        point_gravity(&sun_pos, &pos, univ::MU_SUN)
    }

    fn init(start_time: &AstroTime, stop_time: &AstroTime, settings: &PropSettings) -> ForceSun {
        let duration_days: f64 = *stop_time - *start_time;
        let duration_secs: f64 = duration_days * 86400.0;
        let ntimes = 1 + (duration_secs / settings.gravity_interp_dt_secs).ceil() as u32;
        ForceSun {
            start_time: *start_time,
            sun_interp_dt_secs: settings.sun_interp_dt_secs,
            pos_gcrf: {
                let dt_days = settings.gravity_interp_dt_secs / SECONDS_PER_DAY;
                (0..ntimes)
                    .into_iter()
                    .map(|x| {
                        let tm: AstroTime = *start_time + x as f64 * dt_days;
                        sun::pos_gcrf(&tm)
                    })
                    .collect()
            },
        }
    }
}
