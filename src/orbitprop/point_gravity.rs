use super::types::*;

pub fn point_gravity(body: &Vec3, sat: &Vec3, mu: f64) -> SimpleState {
    let pdiff = body - sat;
    let pdiffnorm = pdiff.norm();
    let bodynorm = body.norm();

    let force = -mu
        * (pdiff / (pdiffnorm * pdiffnorm * pdiffnorm) + body / (bodynorm * bodynorm * bodynorm));

    let mut statedot = SimpleState::zeros();
    statedot
        .index_mut((3.., ..))
        .as_mut_slice()
        .copy_from_slice(force.as_slice());
    statedot
}
