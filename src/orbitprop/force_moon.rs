use super::forceterm::*;
use super::point_gravity::point_gravity;
use super::settings::PropSettings;
use super::types::*;
use super::utils::linterp_idx;

use crate::AstroTime;
use nalgebra as na;

use crate::lpephem::moon;
use crate::univ;

pub struct ForceMoon {
    start_time: AstroTime,
    moon_interp_dt_secs: f64,
    pos_gcrf: Vec<na::Vector3<f64>>,
}

impl ForceTerm<SimpleState> for ForceMoon {
    fn ydot(&self, time: &AstroTime, state: &SimpleState) -> SimpleState {
        let float_idx: f64 = (*time - self.start_time) / self.moon_interp_dt_secs;
        let moon_pos = linterp_idx(&self.pos_gcrf, float_idx).unwrap();
        let pos: Vec3 = Vec3::from_row_slice(state.index((0..3, ..)).as_slice());

        point_gravity(&moon_pos, &pos, univ::MU_MOON)
    }

    fn init(start_time: &AstroTime, stop_time: &AstroTime, settings: &PropSettings) -> ForceMoon {
        let duration_days: f64 = *stop_time - *start_time;
        let duration_secs: f64 = duration_days * 86400.0;
        let ntimes = 1 + (duration_secs / settings.gravity_interp_dt_secs).ceil() as u32;
        ForceMoon {
            start_time: *start_time,
            moon_interp_dt_secs: settings.moon_interp_dt_secs,
            pos_gcrf: {
                (0..ntimes)
                    .into_iter()
                    .map(|x| {
                        let tm: AstroTime = *start_time + x as f64 / 86400.0;
                        moon::pos_gcrf(&tm)
                    })
                    .collect()
            },
        }
    }
}
