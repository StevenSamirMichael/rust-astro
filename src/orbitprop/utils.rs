///
/// Linearly interpolate into a 2-dimensional array
/// given a non-integer floating-point index
///
pub fn linterp_idx<V>(yarr: &Vec<V>, fidx: f64) -> Option<V>
where
    V: std::ops::Mul<f64, Output = V> + std::ops::Add<Output = V> + Copy,
{
    let iidx = fidx as usize;
    if iidx > (yarr.len() - 2) {
        return None;
    }
    let weight: f64 = fidx - (iidx as f64);
    Some((yarr[iidx] * (1.0 - weight)) + (yarr[iidx + 1] * weight))
}
