use super::super::utils::datadir;
use std::io::{self, BufRead};
use std::path::PathBuf;

use super::errors::GravityError;

use nalgebra as na;
type CoeffTable = na::DMatrix<f64>;

type DivisorTable = na::SMatrix<f64, 20, 20>;

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

    // Equations 8-25 through 8-27, Vallado, page 550.
    pub fn accel_vallado<const N: usize>(&self, pos: &Vec3) -> Vec3 {
        let r: f64 = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        let phi: f64 = f64::asin(pos[2] / r);
        let lambda: f64 = f64::atan2(pos[1], pos[0]);
        let p = Gravity::associated_legendre::<N>(phi);
        let tanphi: f64 = f64::tan(phi);

        let reoverr = self.radius / r;

        let mut dudr: f64 = 0.0;
        let mut dudphi: f64 = 0.0;
        let mut dudlambda: f64 = 0.0;

        for l in 0..(N - 1) {
            let reoverr_l = f64::powf(reoverr, l as f64);
            let mut dudr_sum = 0.0;
            let mut dudphi_sum = 0.0;
            let mut dudlambda_sum = 0.0;

            for m in 2..(l + 1) {
                let cosmtheta = f64::cos(m as f64 * lambda);
                let sinmtheta = f64::sin(m as f64 * lambda);
                let clm = self.coeffs[(l, m)];
                let mut slm = 0.0;
                if m > 0 {
                    slm = self.coeffs[(m - 1, l)];
                }
                let m1 = clm * cosmtheta + slm * sinmtheta;
                let m2 = slm * cosmtheta - clm * sinmtheta;

                dudr_sum += p[(l, m)] * m1;
                dudphi_sum += (p[(l, m + 1)] - m as f64 * tanphi * p[(l, m)]) * m1;
                dudlambda_sum += m as f64 * p[(l, m)] * m2;
            }
            dudr_sum *= reoverr_l * (l + 1) as f64;
            dudphi_sum *= reoverr_l;
            dudlambda_sum *= reoverr_l;

            dudr += dudr_sum;
            dudphi += dudphi_sum;
            dudlambda += dudlambda_sum;
        }
        dudr *= -1.0 * self.gravity_constant / r / r;
        dudphi *= self.gravity_constant / r;
        dudlambda *= self.gravity_constant / r;

        //return Vec3::from([dudr, dudphi / r, dudlambda / r]);

        let rho = f64::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
        let term1 = dudr / r - pos[2] / r / r / rho * dudphi;
        let term2 = dudlambda / rho / rho;
        let a1 = term1 * pos[0] - term2 * pos[1];
        let a2 = term1 * pos[1] + term2 * pos[0];
        let a3 = dudr / r * pos[2] + rho / r / r * dudphi;

        Vec3::from([a1, a2, a3]) - pos * self.gravity_constant / r / r / r
    }

    fn associated_legendre<const N: usize>(v: f64) -> na::SMatrix<f64, N, N> {
        let mut p = na::SMatrix::<f64, N, N>::zeros();
        let sinv = v.sin();
        let cosv = v.cos();

        // Set initial values
        p[(0, 0)] = 1.0;
        p[(1, 0)] = sinv;
        p[(1, 1)] = cosv;

        // Walk through to get subsequent values
        // Vallado equation 8-57
        for l in 2..N {
            for m in 0..(l + 1) {
                if m == 0 {
                    p[(l, m)] = ((2 * l - 1) as f64 * p[(l - 1, m)]
                        - (l - 1) as f64 * p[(l - 2, 0)])
                        / (l as f64);
                } else if m < l {
                    p[(l, m)] = p[(l - 2, m)] + (2 * l - 1) as f64 * cosv * p[(l - 1, m - 1)];
                } else {
                    p[(l, l)] = (2 * l - 1) as f64 * cosv * p[(l - 1, l - 1)];
                }
            }
        }

        p
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

    pub fn from_file(filename: &str) -> Result<Gravity, Box<dyn std::error::Error + Send + Sync>> {
        let path = datadir::get().unwrap_or(PathBuf::from(".")).join(filename);
        if !path.is_file() {
            return Err(GravityError::new("File does not exist").into());
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
            return Err(GravityError::new("Invalid File; did not find max_degree").into());
        }

        // Create matrix with lookup values
        let mut cs: CoeffTable = CoeffTable::zeros(max_degree + 1, max_degree + 1);

        for line in &lines[header_cnt..] {
            let s: Vec<&str> = line.split_whitespace().collect();
            if s.len() < 3 {
                println!("bad line = {}", line);
                return Err(GravityError::new("Invalid File").into());
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

    use super::super::super::itrfcoord::ITRFCoord;

    #[test]
    fn load_gravity() {
        use std::f64::consts::PI;
        let g = Gravity::from_file("jgm3.gfc").unwrap();
        let itrf = ITRFCoord::from_geodetic_deg(42.466, -71.1516, 0.0);
        let accel = g.accel_t::<4, 8>(&itrf.into());
        //let truth = [2.34440183, 6.86790166, -6.1888031];
        let accel_ned = itrf.q_ned2itrf().conjugate() * accel;
        println!("ned accel = {}", accel_ned);
        let ew_deflection = f64::atan2(accel_ned[1], accel_ned[2]) * 180.0 / PI * 3600.0;
        let nw_deflection = f64::atan2(accel_ned[0], accel_ned[2]) * 180.0 / PI * 3600.0;
        println!("Deflections = {} :: {}", ew_deflection, nw_deflection);
        let a2 = g.accel_vallado::<8>(&itrf.into());
        println!("accel = {}", accel);
        println!("annorm = {}", accel.norm());

        //println!("a2 = {}", a2);
        //let a2_ned = itrf.q_ned2itrf().conjugate() * a2;
        println!("a2 ned = {}", a2);
        println!("itrf = {}", itrf.itrf);
    }
}
