use crate::utils;
use nalgebra as na;
use std::error::Error;
use std::io::{self, BufRead};
use std::path::PathBuf;

#[derive(Debug)]
pub struct IERSTable {
    data: [Option<na::DMatrix<f64>>; 6],
}

impl IERSTable {
    pub fn from_file(fname: &str) -> Result<IERSTable, Box<dyn Error + Send + Sync>> {
        let mut table = IERSTable {
            data: [None, None, None, None, None, None],
        };

        let path = utils::datadir::get()
            .unwrap_or(PathBuf::from("."))
            .join(fname);
        if !path.is_file() {
            return Err(utils::AstroErr::new("Could not open file").into());
        }

        let mut tnum: i32 = -1;
        let file = std::fs::File::open(&path)?;
        let lines = io::BufReader::new(file).lines();

        for line in lines {
            match line {
                Ok(l) => {
                    let tline = l.trim();
                    if tline.len() < 10 {
                        continue;
                    }
                    if tline[..4].eq("j =") {
                        tnum = tline[5..6].parse()?;
                        if tnum < 0 || tnum > 5 {
                            return Err(utils::AstroErr::new("Error parsing file").into());
                        }
                        table.data[tnum as usize] = Some(na::DMatrix::<f64>::zeros(100, 100));
                        continue;
                    }
                    (&table.data[0].unwrap())[(0, 0)] = 1.0;
                }
                Err(_) => continue,
            }
        }
        Ok(table)
    }
}
