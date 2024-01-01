use nalgebra as na;

use crate::orbitprop;
use crate::orbitprop::PropSettings;
use crate::AstroResult;
use crate::AstroTime;

type PVCovType = na::SMatrix<f64, 6, 6>;

#[derive(Clone, Debug)]
pub enum StateCov {
    None,
    PVCov(PVCovType),
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

    /// Set position uncertainty (1-sigma) in the position, velocity, angular momentum
    /// frame
    pub fn set_pvh_pos_uncertainty(&mut self, sigma_pvh: &na::Vector3<f64>) {
        self.cov = StateCov::PVCov({
            let mut m = PVCovType::zeros();
            let mut diag = na::Vector6::<f64>::zeros();
            diag[0] = sigma_pvh[0] * sigma_pvh[0];
            diag[1] = sigma_pvh[1] * sigma_pvh[1];
            diag[2] = sigma_pvh[2] * sigma_pvh[2];
            m.set_diagonal(&diag);
            m
        })
    }

    pub fn propagate(
        &self,
        time: &AstroTime,
        option_settings: Option<&PropSettings>,
    ) -> AstroResult<SatState> {
        let default = orbitprop::PropSettings::default();
        let settings = option_settings.unwrap_or(&default);
        match self.cov {
            StateCov::None => {
                let res = orbitprop::propagate(&self.pv, &self.time, time, None, settings, None)?;
                Ok(SatState {
                    time: time.clone(),
                    pv: res.state[0],
                    cov: StateCov::None,
                })
            }
            StateCov::PVCov(cov) => {
                let mut state = na::SMatrix::<f64, 6, 7>::zeros();
                state.fixed_view_mut::<6, 1>(0, 0).copy_from(&self.pv);
                let res = orbitprop::propagate(&state, &self.time, time, None, settings, None)?;
                Ok(SatState {
                    time: time.clone(),
                    pv: res.state[0].fixed_view::<6, 1>(0, 0).into(),
                    cov: {
                        let phi = res.state[0].fixed_view::<6, 6>(0, 1);
                        StateCov::PVCov(phi * cov * phi.transpose())
                    },
                })
            }
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            r#"Satellite State
                Time: {}
            Position: [{:+8.0}, {:+8.0}, {:+8.0}] m,
            Velocity: [{:+8.3}, {:+8.3}, {:+8.3}] m/s"#,
            self.time, self.pv[0], self.pv[1], self.pv[2], self.pv[3], self.pv[4], self.pv[5],
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
