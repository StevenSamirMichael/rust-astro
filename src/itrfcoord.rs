use std::f64::consts::PI;
const DEG2RAD: f64 = PI / 180.;
const RAD2DEG: f64 = 180. / PI;

pub const WGS84_A: f64 = 6378137.0;
pub const WGS84_F: f64 = 0.003352810664747;

use nalgebra as na;

use crate::skerror;
use crate::types::Quaternion as Quat;
use crate::types::Vec3;
use crate::SKResult;

///
/// Representation of a coordinate in the
/// International Terrestrial Reference Frame (ITRF)
///
/// This coordinate object can be created from and also
/// output to Geodetic coordinates (latitude, longitude,
/// height above ellipsoid).
///
/// Functions are also available to provide rotation
/// quaternions to the East-North-Up frame
/// and North-East-Down frame at this coordinate
///
#[derive(PartialEq, PartialOrd, Copy, Clone)]
pub struct ITRFCoord {
    pub itrf: Vec3,
}

impl std::fmt::Display for ITRFCoord {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let (lat, lon, hae) = self.to_geodetic_deg();
        write!(
            f,
            "ITRFCoord(lat: {:8.4} deg, lon: {:8.4} deg, altitude: {:5.2} km)",
            lat,
            lon,
            hae / 1.0e3
        )
    }
}

impl std::ops::Add<Vec3> for ITRFCoord {
    type Output = Self;
    fn add(self, other: Vec3) -> Self::Output {
        Self {
            itrf: self.itrf + other,
        }
    }
}

impl std::ops::Add<Vec3> for &ITRFCoord {
    type Output = ITRFCoord;
    fn add(self, other: Vec3) -> Self::Output {
        ITRFCoord {
            itrf: self.itrf + other,
        }
    }
}

impl std::ops::Add<&Vec3> for ITRFCoord {
    type Output = Self;
    fn add(self, other: &Vec3) -> Self::Output {
        Self {
            itrf: self.itrf + other,
        }
    }
}

impl std::ops::Add<&Vec3> for &ITRFCoord {
    type Output = ITRFCoord;
    fn add(self, other: &Vec3) -> Self::Output {
        ITRFCoord {
            itrf: self.itrf + other,
        }
    }
}

impl std::ops::Sub<Vec3> for ITRFCoord {
    type Output = Self;
    fn sub(self, other: Vec3) -> Self::Output {
        Self {
            itrf: self.itrf - other,
        }
    }
}

impl std::ops::Sub<ITRFCoord> for ITRFCoord {
    type Output = Vec3;
    fn sub(self, other: ITRFCoord) -> Vec3 {
        self.itrf - other.itrf
    }
}

impl std::ops::Sub<ITRFCoord> for &ITRFCoord {
    type Output = Vec3;
    fn sub(self, other: ITRFCoord) -> Vec3 {
        self.itrf - other.itrf
    }
}

impl std::ops::Sub<&ITRFCoord> for &ITRFCoord {
    type Output = Vec3;
    fn sub(self, other: &ITRFCoord) -> Vec3 {
        self.itrf - other.itrf
    }
}

impl std::ops::Sub<&ITRFCoord> for ITRFCoord {
    type Output = Vec3;
    fn sub(self, other: &ITRFCoord) -> Vec3 {
        self.itrf - other.itrf
    }
}

impl std::convert::From<[f64; 3]> for ITRFCoord {
    fn from(v: [f64; 3]) -> Self {
        ITRFCoord {
            itrf: Vec3::from(v),
        }
    }
}

impl std::convert::From<&[f64]> for ITRFCoord {
    fn from(v: &[f64]) -> Self {
        assert!(v.len() == 3);
        ITRFCoord {
            itrf: Vec3::from_row_slice(v),
        }
    }
}

impl std::convert::From<Vec3> for ITRFCoord {
    fn from(v: Vec3) -> Self {
        ITRFCoord { itrf: v }
    }
}

impl std::convert::From<ITRFCoord> for Vec3 {
    fn from(itrf: ITRFCoord) -> Self {
        itrf.itrf
    }
}

impl ITRFCoord {
    /// Returns an ITRF Coordinate given the geodetic inputs
    ///   with degree units for latitude & longitude
    ///
    /// # Arguments:
    ///
    /// * `lat` - Geodetic latitude in degrees
    /// * `lon` - Geodetic longitude in degrees
    /// * `hae` - Height above ellipsoid, in meters
    ///
    /// # Examples:
    /// ```
    /// // Create coord for ~ Boston, MA
    /// use astro::itrfcoord::ITRFCoord;
    /// let itrf = ITRFCoord::from_geodetic_deg(42.466, -71.1516, 150.0);
    /// ```
    ///
    pub fn from_geodetic_deg(lat: f64, lon: f64, hae: f64) -> ITRFCoord {
        ITRFCoord::from_geodetic_rad(lat * DEG2RAD, lon * DEG2RAD, hae)
    }

    pub fn from_vector(v: &na::Vector3<f64>) -> ITRFCoord {
        ITRFCoord { itrf: v.clone() }
    }

    /// Return an ITRF coordinate given input slice
    /// representing Cartesian coordinates, in meters
    pub fn from_slice(v: &[f64]) -> SKResult<ITRFCoord> {
        if v.len() != 3 {
            return skerror!("Input slice must have 3 elements");
        }
        Ok(ITRFCoord {
            itrf: Vec3::from_row_slice(v),
        })
    }

    /// Returns an ITRF Coordinate given the geodetic inputs
    ///   with radian units for latitude & longitude
    ///
    /// # Arguments:
    ///
    /// * `lat` - Geodetic latitude in radians
    /// * `lon` - Geodetic longitude in radians
    /// * `hae` - Height above ellipsoid, in meters
    ///
    /// # Examples:
    /// ```
    /// // Create coord for ~ Boston, MA
    /// use astro::itrfcoord::ITRFCoord;
    /// use std::f64::consts::PI;
    /// const DEG2RAD: f64 = PI / 180.0;
    /// let itrf = ITRFCoord::from_geodetic_rad(42.466*DEG2RAD, -71.1516*DEG2RAD, 150.0);
    /// ```
    ///
    pub fn from_geodetic_rad(lat: f64, lon: f64, hae: f64) -> ITRFCoord {
        let sinp: f64 = lat.sin();
        let cosp: f64 = lat.cos();
        let sinl: f64 = lon.sin();
        let cosl: f64 = lon.cos();

        let f2 = (1.0 - WGS84_F).powf(2.0);
        let c = 1.0 / (cosp * cosp + f2 * sinp * sinp).sqrt();
        let s = f2 * c;

        ITRFCoord {
            itrf: Vec3::from([
                (WGS84_A * c + hae) * cosp * cosl,
                (WGS84_A * c + hae) * cosp * sinl,
                (WGS84_A * s + hae) * sinp,
            ]),
        }
    }

    /// Returns 3-element tuple representing geodetic coordinates
    ///
    /// # Tuple contents:
    ///
    /// * `.0` - latitude in radians
    /// * `.1` - longitude in radians
    /// * `.2` - height above ellipsoid, in meters
    ///
    pub fn to_geodetic_rad(&self) -> (f64, f64, f64) {
        const B: f64 = WGS84_A * (1.0 - WGS84_F);
        const E2: f64 = 1.0 - (1.0 - WGS84_F) * (1.0 - WGS84_F);
        const EP2: f64 = E2 / (1.0 - E2);

        let rho = (self.itrf[0] * self.itrf[0] + self.itrf[1] * self.itrf[1]).sqrt();
        let mut beta: f64 = f64::atan2(self.itrf[2], (1.0 - WGS84_F) * rho);
        let mut sinbeta: f64 = beta.sin();
        let mut cosbeta: f64 = beta.cos();
        let mut phi: f64 = f64::atan2(
            self.itrf[2] + B * EP2 * sinbeta.powf(3.0),
            rho - WGS84_A * E2 * cosbeta.powf(3.0),
        );
        let mut betanew: f64 = f64::atan2((1.0 - WGS84_F) * phi.sin(), phi.cos());
        for _x in 0..5 {
            beta = betanew;
            sinbeta = beta.sin();
            cosbeta = beta.cos();
            phi = f64::atan2(
                self.itrf[2] + B * EP2 * sinbeta.powf(3.0),
                rho - WGS84_A * E2 * cosbeta.powf(3.0),
            );
            betanew = f64::atan2((1.0 - WGS84_F) * phi.sin(), phi.cos());
        }
        let lat: f64 = phi;
        let lon: f64 = f64::atan2(self.itrf[1], self.itrf[0]);
        let sinphi: f64 = phi.sin();
        let n: f64 = WGS84_A / (1.0 - E2 * sinphi * sinphi).sqrt();
        let h = rho * phi.cos() + (self.itrf[2] + E2 * n * sinphi) * sinphi - n;
        (lat, lon, h)
    }

    /// Returns 3-element tuple representing geodetic coordinates
    ///
    /// # Tuple contents:
    ///
    /// * `.0` - latitude in degrees
    /// * `.1` - longitude in degrees
    /// * `.2` - height above ellipsoid, in meters
    ///
    pub fn to_geodetic_deg(&self) -> (f64, f64, f64) {
        let (lat_rad, lon_rad, hae) = self.to_geodetic_rad();
        (lat_rad * RAD2DEG, lon_rad * RAD2DEG, hae)
    }

    /// Return geodetic longitude in radians, [-pi, pi]
    ///
    #[inline]
    pub fn longitude_rad(&self) -> f64 {
        f64::atan2(self.itrf[1], self.itrf[0])
    }

    /// Return geodetic longitude in degrees, [-180, 180]
    #[inline]
    pub fn longitude_deg(&self) -> f64 {
        self.longitude_rad() * RAD2DEG
    }

    /// return geodetic latitude in radians, [-pi/2, pi/2]
    #[inline]
    pub fn latitude_rad(&self) -> f64 {
        let (lat, _a, _b) = self.to_geodetic_rad();
        lat
    }

    /// Return height above ellipsoid in meters
    #[inline]
    pub fn hae(&self) -> f64 {
        let (_a, _b, hae) = self.to_geodetic_rad();
        hae
    }

    /// Return geodetic latitude in degrees, [-pi/2, pi/2]
    #[inline]
    pub fn latitude_deg(&self) -> f64 {
        self.latitude_rad() * RAD2DEG
    }

    /// Return quaternion representing rotation from the
    /// North-East-Down (NED) coordinate frame to the
    /// ITRF coordinate frame
    #[inline]
    pub fn q_ned2itrf(&self) -> Quat {
        let (lat, lon, _) = self.to_geodetic_rad();
        Quat::from_axis_angle(&Vec3::z_axis(), lon)
            * Quat::from_axis_angle(&Vec3::y_axis(), -lat - PI / 2.0)
    }

    /// Convert coordinate to a North-East-Down (NED)
    /// coordinate relative to a reference coordinate
    ///
    /// # Arguemnts
    ///
    /// * ref_coord - &ITRFCoord representing reference
    ///
    /// # Return
    ///
    /// * nalgebra::Vector3<f64> representing NED position
    ///   relative to reference.  Units are meters
    ///
    /// # Examples:
    /// ```
    /// use astro::itrfcoord::ITRFCoord;
    /// // Create coord
    /// let itrf1 = ITRFCoord::from_geodetic_deg(42.466, -71.1516, 150.0);
    /// // Crate 2nd coord 100 meters above
    /// let itrf2 = ITRFCoord::from_geodetic_deg(42.466, -71.1516, 250.0);
    ///
    /// // Get NED of itrf1 relative to itrf2
    /// let ned = itrf1.to_ned(&itrf2);
    /// // Should return [0.0, 0.0, 100.0]
    /// ```
    ///
    pub fn to_ned(&self, ref_coord: &ITRFCoord) -> Vec3 {
        self.q_ned2itrf().conjugate() * (self.itrf - ref_coord.itrf)
    }

    /// Return quaternion representing rotation from the
    /// East-North-Up (ENU) coordinate frame to the
    /// ITRF coordinate frame
    pub fn q_enu2itrf(&self) -> Quat {
        let (lat, lon, _) = self.to_geodetic_rad();
        Quat::from_axis_angle(&Vec3::z_axis(), lon + PI / 2.0)
            * Quat::from_axis_angle(&Vec3::x_axis(), PI / 2.0 - lat)
    }

    /// Convert coordinate to a East-North-Up (ENU)
    /// coordinate relative to a reference coordinate
    ///
    /// # Arguemnts
    ///
    /// * ref_coord - &ITRFCoord representing reference
    ///
    /// # Return
    ///
    /// * nalgebra::Vector3<f64> representing ENU position
    ///   relative to reference.  Units are meters
    ///
    /// # Examples:
    /// ```
    /// use astro::itrfcoord::ITRFCoord;
    /// // Create coord
    /// let itrf1 = ITRFCoord::from_geodetic_deg(42.466, -71.1516, 150.0);
    /// // Crate 2nd coord 100 meters above
    /// let itrf2 = ITRFCoord::from_geodetic_deg(42.466, -71.1516, 250.0);
    ///
    /// // Get ENU of itrf1 relative to itrf2
    /// let enu = itrf1.to_ned(&itrf2);
    /// // Should return [0.0, 0.0, -100.0]
    /// ```
    ///
    pub fn to_enu(&self, other: &ITRFCoord) -> Vec3 {
        self.q_enu2itrf().conjugate() * (self.itrf - other.itrf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geodetic() {
        let lat_deg: f64 = 42.466;
        let lon_deg: f64 = -71.0;
        let hae: f64 = 150.0;
        let itrf = ITRFCoord::from_geodetic_deg(lat_deg, lon_deg, hae);
        println!("{}", itrf);
        // Check conversions
        assert!(((lat_deg - 42.466) / 42.466).abs() < 1.0e-6);
        assert!(((lon_deg + 71.0) / 71.0).abs() < 1.0e-6);
        assert!(((hae - 150.0) / 150.0).abs() < 1.0e-6);
    }

    #[test]
    fn test_ned_enu() {
        let lat_deg: f64 = 42.466;
        let lon_deg: f64 = -74.0;
        let hae: f64 = 150.0;
        let itrf1 = ITRFCoord::from_geodetic_deg(lat_deg, lon_deg, hae);
        let itrf2 = ITRFCoord::from_geodetic_deg(lat_deg, lon_deg, hae + 100.0);
        let ned = itrf2.to_ned(&itrf1);
        let enu = itrf2.to_enu(&itrf1);
        assert!(enu[0].abs() < 1.0e-6);
        assert!(enu[1].abs() < 1.0e-6);
        assert!(((enu[2] - 100.0) / 100.0).abs() < 1.0e-6);
        assert!(ned[0].abs() < 1.0e-6);
        assert!(ned[1].abs() < 1.0e-6);
        assert!(((ned[2] + 100.0) / 100.0).abs() < 1.0e-6);

        let dvec = Vec3::from([-100.0, -200.0, 300.0]);
        let itrf3 = itrf2 + itrf2.q_ned2itrf() * dvec;
        let nedvec = itrf3.to_ned(&itrf2);
        let itrf4 = itrf2 + itrf2.q_enu2itrf() * dvec;
        let enuvec = itrf4.to_enu(&itrf2);
        for x in 0..3 {
            assert!(((nedvec[x] - dvec[x]) / nedvec[x]).abs() < 1.0e-3);
            assert!(((enuvec[x] - dvec[x]) / nedvec[x]).abs() < 1.0e-3);
        }
        /*
        let q = Quat::from_axis_angle(&Vec3::z_axis(), -0.003);
        println!("{}", q);
        println!("{}", q.to_rotation_matrix());
        */

        let itrf1 = ITRFCoord::from_geodetic_deg(lat_deg, lon_deg, hae);
        let itrf2 = itrf1 + itrf1.q_ned2itrf() * na::vector![0.0, 0.0, 10000.0];
        println!("height diff = {}", itrf2.hae() - itrf1.hae());
    }
}
