use nalgebra as na;

use crate::orbitprop;
use crate::orbitprop::PropSettings;
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

    pub fn propagate(
        &self,
        time: &AstroTime,
        option_settings: Option<&PropSettings>,
    ) -> AstroResult<SatState> {
        let default = orbitprop::PropSettings::default();
        let settings = option_settings.unwrap_or(&default);
        let res = orbitprop::propagate(&self.pv, &self.time, time, None, settings, None)?;
        Ok(SatState {
            time: time.clone(),
            pv: res.state[0],
            cov: StateCov::None,
        })
    }

    pub fn to_string(&self) -> String {
        format!(
            r#"Satellite State
                Time: {}
            Position: {}
            Velocity: {}"#,
            self.time,
            self.pv.fixed_view::<3, 1>(0, 0).transpose(),
            self.pv.fixed_view::<3, 1>(3, 0).transpose(),
        )
    }
}

impl std::fmt::Display for SatState {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
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

        let state2 = satstate.propagate(&(satstate.time + 1.0), None)?;

        println!("state 2 = {:?}", state2);

        let state0 = state2.propagate(&satstate.time, None);
        println!("state 0 = {:?}", state0);
        Ok(())
    }
}
