use super::super::datadir;
use std::io::{self, BufRead};
use std::path::PathBuf;

use super::errors::GravityError;

use nalgebra as na;
type CoeffTable = na::DMatrix<f64>;

const MAX_ORDER: usize = 16;
type DivisorTable = na::SMatrix<f64, 20, 20>;

pub struct Gravity {
    pub name: String,
    pub gravity_constant: f64,
    pub radius: f64,
    pub max_degree: usize,
    pub coeffs: CoeffTable,
    pub divisor_table: DivisorTable,
}

impl Gravity {
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
            header_cnt = header_cnt + 1;
        }
        if max_degree == 0 {
            return Err(GravityError::new("Invalid File").into());
        }

        // Create matrix with lookup values
        let mut cs: CoeffTable = CoeffTable::zeros(max_degree + 1, max_degree + 1);

        for line in &lines[header_cnt..] {
            let s: Vec<&str> = line.split_whitespace().collect();
            if s.len() < 3 {
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
                let mut scale = 1.0f64;
                for k in (n - m + 1)..(n + m + 1) {
                    scale = scale * k as f64;
                }
                scale = scale / (2.0 * n as f64 + 1.0);
                scale = f64::sqrt(scale);
                cs[(n, m)] = cs[(n, m)] * scale;

                if m > 0 {
                    cs[(m - 1, n)] = cs[(m - 1, n)] * scale;
                }
            }
        }

        for m in 0..19 {}

        Ok(Gravity {
            name: String::from(name),
            gravity_constant: gravity_constant,
            radius: radius,
            max_degree: max_degree,
            coeffs: cs,
            divisor_table: {
                let mut dt: DivisorTable = DivisorTable::zeros();
                for m in 0..19 {
                    dt[(m, m)] = 2.0 * m as f64 - 1.0;
                    let n = m + 1;
                    dt[(n, m)] = (2.0 * n as f64 - 1.0) / (n - m) as f64;
                    for n in (m + 2)..19 {
                        dt[(n, m)] = (2.0 * n as f64 - 1.0) / (n - m) as f64;
                    }
                }
                dt
            },
        })
    }
}
