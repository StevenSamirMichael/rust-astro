pub use crate::frametransform::Quat;
pub use crate::utils::AstroResult;

use nalgebra as na;
pub type Vec3 = na::Vector3<f64>;

pub type SimpleState = na::Vector6<f64>;
pub type CovState = na::Matrix<f64, na::Const<6>, na::Const<7>, na::ArrayStorage<f64, 6, 7>>;
