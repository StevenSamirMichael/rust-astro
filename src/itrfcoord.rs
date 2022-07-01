extern crate nalgebra;

use std::f64::consts::PI;
const DEG2RAD: f64 = PI / 180.;
const RAD2DEG: f64 = 180. / PI;

const WGS84_A: f64 = 6378137.0;
const WGS84_F: f64 = 0.003352810664747;

use nalgebra as na;
type Vec3 = na::Vector3<f64>;
type Quat = na::UnitQuaternion<f64>;

#[derive(PartialEq, PartialOrd, Copy, Clone)]
pub struct ITRFCoord {
    itrf: Vec3,
}

impl std::fmt::Display for ITRFCoord {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let (lat, lon, hae) = self.to_geodetic_deg();
        write!(
            f,
            "ITRFCoord(lat: {:8.4} deg, lon: {:8.4} deg, hae: {:5.2} km)",
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
    pub fn from_geodetic_deg(lat: f64, lon: f64, hae: f64) -> ITRFCoord {
        ITRFCoord::from_geodetic_rad(lat * DEG2RAD, lon * DEG2RAD, hae)
    }

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

    pub fn to_geodetic_deg(&self) -> (f64, f64, f64) {
        let (lat_rad, lon_rad, hae) = self.to_geodetic_rad();
        (lat_rad * RAD2DEG, lon_rad * RAD2DEG, hae)
    }

    #[inline]
    pub fn longitude_rad(&self) -> f64 {
        f64::atan2(self.itrf[1], self.itrf[0])
    }

    #[inline]
    pub fn longitude_deg(&self) -> f64 {
        self.longitude_rad() * RAD2DEG
    }

    #[inline]
    pub fn latitude_rad(&self) -> f64 {
        let (lat, _a, _b) = self.to_geodetic_rad();
        lat
    }

    #[inline]
    pub fn hae(&self) -> f64 {
        let (_a, _b, hae) = self.to_geodetic_rad();
        hae
    }

    #[inline]
    pub fn latitude_deg(&self) -> f64 {
        self.latitude_rad() * RAD2DEG
    }

    #[inline]
    pub fn q_ned2itrf(&self) -> Quat {
        let (lat, lon, _) = self.to_geodetic_rad();
        Quat::from_axis_angle(&Vec3::z_axis(), lon)
            * Quat::from_axis_angle(&Vec3::y_axis(), -lat - PI / 2.0)
    }

    pub fn to_ned(&self, other: &ITRFCoord) -> Vec3 {
        self.q_ned2itrf().conjugate() * (self.itrf - other.itrf)
    }

    pub fn q_enu2itrf(&self) -> Quat {
        let (lat, lon, _) = self.to_geodetic_rad();
        Quat::from_axis_angle(&Vec3::z_axis(), lon + PI / 2.0)
            * Quat::from_axis_angle(&Vec3::x_axis(), PI / 2.0 - lat)
    }

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
    }
}
