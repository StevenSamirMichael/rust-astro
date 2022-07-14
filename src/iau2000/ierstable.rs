use crate::utils::{self, AstroErr, AstroResult};
use nalgebra as na;
use std::io::{self, BufRead};
use std::path::PathBuf;

#[derive(Debug)]
pub struct IERSTable {
    data: [na::DMatrix<f64>; 6],
}

impl IERSTable {
    pub fn from_file(fname: &str) -> AstroResult<IERSTable> {
        let mut table = IERSTable {
            data: [
                na::DMatrix::<f64>::zeros(0, 0),
                na::DMatrix::<f64>::zeros(0, 0),
                na::DMatrix::<f64>::zeros(0, 0),
                na::DMatrix::<f64>::zeros(0, 0),
                na::DMatrix::<f64>::zeros(0, 0),
                na::DMatrix::<f64>::zeros(0, 0),
            ],
        };

        let path = utils::datadir::get()
            .unwrap_or(PathBuf::from("."))
            .join(fname);
        if !path.is_file() {
            return utils::astroerr!("Could not open file: {}", fname);
        }

        let mut tnum: i32 = -1;
        let mut rowcnt: usize = 0;
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
                        let s: Vec<&str> = tline.split_whitespace().collect();
                        let tsize: usize = s[s.len() - 1].parse().unwrap_or(0);
                        if tnum < 0 || tnum > 5 || tsize == 0 {
                            return utils::astroerr!(
                                "Error parsing file {}, invalid table definition line",
                                fname
                            );
                        }
                        table.data[tnum as usize - 1] = na::DMatrix::<f64>::zeros(tsize, 17);
                        rowcnt = 0;
                        continue;
                    } else if tnum >= 0 {
                        if table.data[tnum as usize - 1].nrows() < 17 {
                            return Err(utils::AstroErr::new(
                                format!("Error parsing file {}, table not initialized", fname)
                                    .as_str(),
                            )
                            .into());
                        }
                        table.data[tnum as usize - 1].set_row(
                            rowcnt,
                            &na::SMatrix::<f64, 1, 17>::from_iterator(
                                tline
                                    .split_whitespace()
                                    .into_iter()
                                    .map(|x| x.parse().unwrap()),
                            ),
                        );
                        rowcnt = rowcnt + 1;
                    }
                }
                Err(_) => continue,
            }
        }
        Ok(table)
    }
}

#[cfg(test)]
mod tests {
    use super::IERSTable;

    #[test]
    fn load_table() {
        let t = IERSTable::from_file("tab5.2a.txt");
        if t.is_err() {
            panic!("{}", t.unwrap_err());
        }
        println!("got t: {}", t.is_ok());
    }
}
