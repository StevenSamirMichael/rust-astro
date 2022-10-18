//pub fn linterp<T, V>(xarr: &Vec<T>, yarr: &Vec<V>, xval: T) -> Option<V>
//where
//    T: PartialOrd + std::ops::Sub<Output = f64> + Copy,
//    V: std::ops::Mul<f64, Output = V> + std::ops::Add<Output = V> + Copy,
//{
//    match xarr.into_iter().position(|x| *x > xval) {
//        None => None,
//        Some(idx) => {
//            if idx == 0 {
//                return None;
//            }
//            if idx > (yarr.len() - 1) {
//                return None;
//            }
//            let weight: f64 = (xarr[idx] - xval) / (xarr[idx] - xarr[idx - 1]);
//            Some((yarr[idx] * (1.0 - weight)) + (yarr[idx - 1] * weight))
//        }
//    }
//}

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
