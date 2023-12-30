use nalgebra as na;

use crate::orbitprop;
use crate::AstroResult;
use crate::AstroTime;

#[derive(Clone, Debug)]
pub enum StateCov {
    PVCov(na::SMatrix<f64, 6, 6>),
    PVRCov(na::SMatrix<f64, 6, 7>),
    None,
}

#[derive(Clone, Debug)]
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

    pub fn propagate(&self, time: &AstroTime) -> AstroResult<SatState> {
        let mut settings = orbitprop::PropSettings::default();
        settings.gravity_order = 2;
        let res = orbitprop::propagate(&self.pv, &self.time, time, None, &settings, None)?;
        Ok(SatState {
            time: time.clone(),
            pv: res.state[0],
            cov: StateCov::None,
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::univ;

    #[test]
    fn test_satstate() -> AstroResult<()> {
        let satstate = SatState::from_pv(
            &AstroTime::from_datetime(2015, 3, 20, 0, 0, 0.0),
            &na::vector![univ::GEO_R, 0.0, 0.0],
            &na::vector![0.0, (univ::MU_EARTH / univ::GEO_R).sqrt(), 0.0],
        );
        println!("state orig = {:?}", satstate);

        let state2 = satstate.propagate(&(satstate.time + 1.0))?;

        println!("state 2 = {:?}", state2);

        let state0 = state2.propagate(&satstate.time);
        println!("state 0 = {:?}", state0);
        Ok(())
    }
}
