use super::forceterm::*;
use super::settings::PropSettings;
use super::types::*;
use super::utils::linterp_idx;
use crate::astrotime::AstroTime;
use crate::frametransform as ft;

use crate::gravity::GRAVITY_JGM3;
use nalgebra as na;

pub struct ForceGravity {
    gravity_order: usize,
    gravity_interp_dt_secs: f64,
    start_time: AstroTime,
    qitrf2gcrf_arr: Vec<na::Quaternion<f64>>,
}

impl ForceTerm<SimpleState> for ForceGravity {
    fn ydot(&self, time: &AstroTime, state: &SimpleState) -> SimpleState {
        let float_idx: f64 = (*time - self.start_time) / self.gravity_interp_dt_secs;
        let qitrf2gcrf: Quat = na::UnitQuaternion::from_quaternion(
            linterp_idx(&self.qitrf2gcrf_arr, float_idx).unwrap(),
        );
        let pos: Vec3 = Vec3::from_row_slice(state.index((0..3, ..)).as_slice());
        let pitrf: Vec3 = qitrf2gcrf.conjugate() * pos;
        let accel: Vec3 = GRAVITY_JGM3.accel(&pitrf, self.gravity_order);
        let accel_gcrf: Vec3 = qitrf2gcrf * accel;

        let mut statedot = SimpleState::zeros();
        statedot
            .index_mut((3.., ..))
            .as_mut_slice()
            .copy_from_slice(accel_gcrf.as_slice());
        statedot
    }

    fn init(
        start_time: &AstroTime,
        stop_time: &AstroTime,
        settings: &PropSettings,
    ) -> ForceGravity {
        let duration_days: f64 = *stop_time - *start_time;
        let duration_secs: f64 = duration_days * 86400.0;
        ForceGravity {
            gravity_order: settings.gravity_order as usize,
            gravity_interp_dt_secs: settings.gravity_interp_dt_secs,
            start_time: *start_time,
            qitrf2gcrf_arr: {
                let ntimes = 1 + (duration_secs / settings.gravity_interp_dt_secs).ceil() as u32;
                let ret = (0..ntimes)
                    .into_iter()
                    .map(|x| {
                        let tm = *start_time + x as f64 / 86400.0;
                        ft::qgcrf2itrf_approx(&tm).quaternion().to_owned()
                    })
                    .collect();
                ret
            },
        }
    }
}
