use nalgebra as na;

use crate::AstroTime;

pub enum StateCov {
    PVCov(na::SMatrix<f64, 6, 6>),
    PVRCov(na::SMatrix<f64, 6, 7>),
    None,
}

pub struct SatState {
    pub time: AstroTime,
    pub pv: na::Vector6<f64>,
    pub cov: StateCov,
}

impl SatState {
    pub fn from_pv(time: &AstroTime, pos: &na::Vector3<f64>, vel: &na::Vector3<f64>) -> SatState {
        SatState {
            time: time.clone(),
            pv: na::vector![pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]],
            cov: StateCov::None,
        }
    }
}
