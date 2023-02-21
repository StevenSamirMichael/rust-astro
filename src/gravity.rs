use crate::utils::*;
use std::collections::HashMap;
use std::io::{self, BufRead};
use std::path::PathBuf;

use nalgebra as na;
type CoeffTable = na::DMatrix<f64>;

type DivisorTable = na::SMatrix<f64, 20, 20>;

///
/// Gravity model enumeration
///
/// For details of models, see:
/// http://icgem.gfz-potsdam.de/tom_longtime
///
#[derive(PartialEq, Eq, Hash)]
pub enum GravityModel {
    JGM3,
    JGM2,
    EGM96,
    ITUGrace16,
}

// Since gravity models don't change, they work well as
// global variables.  Declare them here; they don't actually
// get instantiated until first use.
lazy_static::lazy_static! {
    pub static ref GRAVITY_JGM3: Gravity = Gravity::from_file("JGM3.gfc").unwrap();
    pub static ref GRAVITY_JGM2: Gravity = Gravity::from_file("JGM2.gfc").unwrap();
    pub static ref GRAVITY_EGM96: Gravity = Gravity::from_file("EGM96.gfc").unwrap();
    pub static ref GRAVITY_ITUGRACE16: Gravity = Gravity::from_file("ITU_GRACE16.gfc").unwrap();

    static ref GRAVHASH: HashMap<GravityModel, &'static Gravity> = {
        let mut m = HashMap::new();
        let g1: &Gravity = &GRAVITY_JGM3;
        m.insert(GravityModel::JGM3, g1);
        m.insert(GravityModel::JGM2, &GRAVITY_JGM2);
        m.insert(GravityModel::EGM96, &GRAVITY_EGM96);
        m.insert(GravityModel::ITUGrace16, &GRAVITY_ITUGRACE16);
        m
    };
}

///
/// Return acceleration due to Earth gravity at the input position. The
/// acceleration does not include the centrifugal force, and is output
/// in m/s^2 in the International Terrestrial Reference Frame (ITRF)
///
/// Inputs:
///
///       pos:    Position as nalgebra 3-vecotr
///
///     order:    The order of the gravity model to use.
///               Maximum is 16
///
///     model:    The gravity model to use, of type "GravityModel"
///
///               For details of models, see:
///               http://icgem.gfz-potsdam.de/tom_longtime
///
///               For details of calculation, see Chapter 3.2 of:
///               "Satellite Orbits: Models, Methods, Applications",
///               O. Montenbruck and B. Gill, Springer, 2012.
///
pub fn accel(pos_itrf: &Vec3, order: usize, model: GravityModel) -> Vec3 {
    GRAVHASH.get(&model).unwrap().accel(pos_itrf, order)
}

pub fn accel_jgm3(pos_itrf: &Vec3, order: usize) -> Vec3 {
    GRAVITY_JGM3.accel(pos_itrf, order)
}

#[derive(Debug, Clone)]
pub struct Gravity {
    pub name: String,
    pub gravity_constant: f64,
    pub radius: f64,
    pub max_degree: usize,
    pub coeffs: CoeffTable,
    pub divisor_table: DivisorTable,
    pub divisor_table2: DivisorTable,
}

type Legendre<const N: usize> = na::SMatrix<f64, N, N>;
type Vec3 = na::Vector3<f64>;

///
/// Return acceleration due to Earth gravity at the input position. The
/// acceleration does not include the centrifugal force, and is output
/// in m/s^2 in the International Terrestrial Reference Frame (ITRF)
///
/// Inputs:
///
///       pos:   Position as ITRF coordinate (astro.itrfcoord) or numpy
///              3-vector representing ITRF position in meters
///
///     order:   Order of the gravity model
///
/// See Equation 3.33 of Montenbruck & Gill (referenced above) for
/// calculation details.
///

impl Gravity {
    pub fn accel(&self, pos: &Vec3, order: usize) -> Vec3 {
        // This is tedious, but using generics allows for vectors to be
        // allocated on the stack, which is faster
        if order == 2 {
            self.accel_t::<2, 6>(pos)
        } else if order == 3 {
            self.accel_t::<3, 7>(pos)
        } else if order == 4 {
            self.accel_t::<4, 8>(pos)
        } else if order == 5 {
            self.accel_t::<5, 9>(pos)
        } else if order == 6 {
            self.accel_t::<6, 10>(pos)
        } else if order == 7 {
            self.accel_t::<7, 11>(pos)
        } else if order == 8 {
            self.accel_t::<8, 12>(pos)
        } else if order == 9 {
            self.accel_t::<9, 13>(pos)
        } else if order == 10 {
            self.accel_t::<10, 14>(pos)
        } else if order == 11 {
            self.accel_t::<11, 15>(pos)
        } else if order == 12 {
            self.accel_t::<12, 16>(pos)
        } else if order == 13 {
            self.accel_t::<13, 17>(pos)
        } else if order == 14 {
            self.accel_t::<14, 18>(pos)
        } else if order == 15 {
            self.accel_t::<15, 19>(pos)
        } else {
            self.accel_t::<16, 20>(pos)
        }
    }

    fn accel_t<const N: usize, const NP4: usize>(&self, pos: &Vec3) -> Vec3 {
        let (v, w) = self.compute_legendre::<NP4>(pos);
        let accel = self.accel_from_legendre_t::<N, NP4>(&v, &w);
        accel
    }

    /// See Equation 3.33 in Montenbruck & Gill
    fn accel_from_legendre_t<const N: usize, const NP4: usize>(
        &self,
        v: &Legendre<NP4>,
        w: &Legendre<NP4>,
    ) -> Vec3 {
        let mut accel = Vec3::zeros();

        for n in 0..(N + 1) {
            for m in 0..(n + 1) {
                let cnm = self.coeffs[(n, m)];
                let mut snm = 0.0;
                if m > 0 {
                    snm = self.coeffs[(m - 1, n)];
                }
                if m == 0 {
                    accel[0] -= cnm * v[(n + 1, 1)];
                    accel[1] -= cnm * w[(n + 1, 1)];
                } else {
                    accel[0] += 0.5
                        * ((-cnm * v[(n + 1, m + 1)] - snm * w[(n + 1, m + 1)])
                            + (n - m + 2) as f64
                                * (n - m + 1) as f64
                                * (cnm * v[(n + 1, m - 1)] + snm * w[(n + 1, m - 1)]));

                    accel[1] += 0.5
                        * (-cnm * w[(n + 1, m + 1)]
                            + snm * v[(n + 1, m + 1)]
                            + (n - m + 2) as f64
                                * (n - m + 1) as f64
                                * (-1.0 * cnm * w[(n + 1, m - 1)] + snm * v[(n + 1, m - 1)]));
                }
                accel[2] += (n - m + 1) as f64 * (-1.0 * cnm * v[(n + 1, m)] - snm * w[(n + 1, m)]);
            }
        }

        accel * self.gravity_constant / self.radius / self.radius
    }

    fn compute_legendre<const NP4: usize>(&self, pos: &Vec3) -> (Legendre<NP4>, Legendre<NP4>) {
        let rsq = pos.norm_squared();
        let xfac = pos[0] * self.radius / rsq;
        let yfac = pos[1] * self.radius / rsq;
        let zfac = pos[2] * self.radius / rsq;
        let rfac = self.radius * self.radius / rsq;

        let mut v = Legendre::<NP4>::zeros();
        let mut w = Legendre::<NP4>::zeros();

        v[(0, 0)] = self.radius / rsq.sqrt();
        w[(0, 0)] = 0.0;

        for m in 0..(NP4 - 1) {
            // Along diagnoals
            if m > 0 {
                v[(m, m)] = self.divisor_table[(m, m)]
                    * (xfac * v[(m - 1, m - 1)] - yfac * w[(m - 1, m - 1)]);
                w[(m, m)] = self.divisor_table[(m, m)]
                    * (xfac * w[(m - 1, m - 1)] + yfac * v[(m - 1, m - 1)]);
            }
            // Work down tree
            let n = m + 1;
            v[(n, m)] = self.divisor_table[(n, m)] * zfac * v[(n - 1, m)];
            w[(n, m)] = self.divisor_table[(n, m)] * zfac * w[(n - 1, m)];

            for n in (m + 2)..(NP4 - 1) {
                v[(n, m)] = self.divisor_table[(n, m)] * zfac * v[(n - 1, m)]
                    - self.divisor_table2[(n, m)] * rfac * v[(n - 2, m)];

                w[(n, m)] = self.divisor_table[(n, m)] * zfac * w[(n - 1, m)]
                    - self.divisor_table2[(n, m)] * rfac * w[(n - 2, m)];
            }
        }
        (v, w)
    }

    /// Load Gravity model coefficients from file
    /// Files are at:
    /// http://icgem.gfz-potsdam.de/tom_longtime
    pub fn from_file(filename: &str) -> AstroResult<Gravity> {
        let path = datadir::get().unwrap_or(PathBuf::from(".")).join(filename);
        if !path.is_file() {
            return astroerr!("File does not exist");
        }
        let file = std::fs::File::open(&path)?;

        let mut name = String::new();
        let mut gravity_constant: f64 = 0.0;
        let mut radius: f64 = 0.0;
        let mut max_degree: usize = 0;
        let mut header_cnt = 0;

        let lines: Vec<String> = io::BufReader::new(file)
            .lines()
            .map(|x| x.unwrap_or(String::from("")))
            .collect();

        // Read header lines
        for line in &lines {
            header_cnt = header_cnt + 1;

            let s: Vec<&str> = line.split_whitespace().collect();
            if s.len() < 2 {
                continue;
            }
            if s[0] == "modelname" {
                name = String::from(s[1]);
            } else if s[0] == "earth_gravity_constant" {
                gravity_constant = s[1].parse::<f64>()?;
            } else if s[0] == "radius" {
                radius = s[1].parse::<f64>()?;
            } else if s[0] == "max_degree" {
                max_degree = s[1].parse::<usize>()?;
                //cs = Some(na::DMatrix::<f64>::zeros(
                //    (max_degree + 1) as usize,
                //    (max_degree + 1) as usize,
                //));
            } else if s[0] == "end_of_head" {
                break;
            }
        }
        if max_degree == 0 {
            return astroerr!("Invalid file; did not find max degree");
        }

        // Create matrix with lookup values
        let mut cs: CoeffTable = CoeffTable::zeros(max_degree + 1, max_degree + 1);

        for line in &lines[header_cnt..] {
            let s: Vec<&str> = line.split_whitespace().collect();
            if s.len() < 3 {
                return astroerr!("Invalid line: {}", line);
            }

            let n: usize = s[1].parse()?;
            let m: usize = s[2].parse()?;
            let v1: f64 = s[3].parse()?;
            cs[(n, m)] = v1;
            if m > 0 {
                let v2: f64 = s[4].parse()?;
                cs[(m - 1, n)] = v2;
            }
        }

        // Convert from normalized coefficients to actual coefficients
        for n in 0..(max_degree + 1) {
            for m in 0..(n + 1) {
                let mut scale: f64 = 1.0;
                for k in (n - m + 1)..(n + m + 1) {
                    scale = scale * k as f64;
                }
                scale = scale / (2.0 * n as f64 + 1.0);
                if m > 0 {
                    scale = scale / 2.0;
                }
                scale = 1.0 / f64::sqrt(scale);
                cs[(n, m)] = cs[(n, m)] * scale;

                if m > 0 {
                    cs[(m - 1, n)] = cs[(m - 1, n)] * scale;
                }
            }
        }

        Ok(Gravity {
            name: String::from(name),
            gravity_constant: gravity_constant,
            radius: radius,
            max_degree: max_degree,
            coeffs: cs,
            divisor_table: {
                let mut dt: DivisorTable = DivisorTable::zeros();
                for m in 0..19 {
                    if m > 0 {
                        dt[(m, m)] = 2.0 * m as f64 - 1.0;
                    }
                    let n = m + 1;
                    dt[(n, m)] = (2.0 * n as f64 - 1.0) / (n - m) as f64;
                    for n in (m + 2)..19 {
                        dt[(n, m)] = (2.0 * n as f64 - 1.0) / (n - m) as f64;
                    }
                }
                dt
            },
            divisor_table2: {
                let mut dt: DivisorTable = DivisorTable::zeros();
                for m in 0..19 {
                    for n in (m + 2)..19 {
                        dt[(n, m)] = (n as f64 + m as f64 - 1.0) / (n - m) as f64;
                    }
                }
                dt
            },
        })
    }
}

#[cfg(test)]

mod tests {

    use super::Gravity;

    use crate::itrfcoord::ITRFCoord;
    use crate::types::Vec3;
    use crate::univ::OMEGA_EARTH;
    use std::f64::consts::PI;

    #[test]
    fn test_gravity() {
        // Lexington, ma
        let latitude: f64 = 42.4473;
        let longitude: f64 = -71.2272;
        let altitude: f64 = 0.0;

        // reference gravity computations, using
        // JGM3 model, with 16 terms, found at:
        // http://icgem.gfz-potsdam.de/calcstat/
        // Outputs from above web page below:
        let reference_gravitation: f64 = 9.822206169031;
        // "gravity" includes centrifugal force, "gravitation" does not
        let reference_gravity: f64 = 9.803696372738;
        // Gravity deflections from normal along east-west and north-south
        // direction, in arcseconds
        let reference_ew_deflection_asec: f64 = -1.283542043355E+00;
        let reference_ns_deflection_asec: f64 = -1.311709802440E+00;

        let g = Gravity::from_file("JGM3.gfc").unwrap();
        let coord = ITRFCoord::from_geodetic_deg(latitude, longitude, altitude);
        let gravitation: Vec3 = g.accel(&coord.into(), 16);
        let centrifugal: Vec3 =
            Vec3::new(coord.itrf[0], coord.itrf[1], 0.0) * OMEGA_EARTH * OMEGA_EARTH;
        let gravity = gravitation + centrifugal;

        // Check gravitation matches the reference value
        // from http://icgem.gfz-potsdam.de/calcstat/
        assert!(f64::abs(gravitation.norm() / reference_gravitation - 1.0) < 1.0E-9);
        // Check that gravity matches reference value
        assert!(f64::abs(gravity.norm() / reference_gravity - 1.0) < 1.0E-9);

        // Rotate to ENU coordinate frame
        let g_enu: Vec3 = coord.q_enu2itrf().conjugate() * gravity;

        // Compute East/West and North/South deflections, in arcsec
        let ew_deflection: f64 = -f64::atan2(g_enu[0], -g_enu[2]) * 180.0 / PI * 3600.0;
        let ns_deflection: f64 = -f64::atan2(g_enu[1], -g_enu[2]) * 180.0 / PI * 3600.0;

        // Compare with reference values
        assert!(f64::abs(ew_deflection / reference_ew_deflection_asec - 1.0) < 1.0e-5);
        assert!(f64::abs(ns_deflection / reference_ns_deflection_asec - 1.0) < 1.0e-5);
    }
}
