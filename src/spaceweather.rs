use std::cmp::Ordering;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;

use crate::astrotime::AstroTime;
use crate::utils::{datadir, download_file, skerror, testdirs, SKResult};

use std::sync::RwLock;

use once_cell::sync::OnceCell;

#[derive(Debug, Clone)]
pub struct SpaceWeatherRecord {
    /// Date of record
    pub date: AstroTime,
    /// Bartels Solar Radiation Number.  
    /// A sequence of 27-day intervals counted continuously from 1832 February 8
    pub bsrn: i32,
    /// Number of day within the bsrn
    pub nd: i32,
    /// Kp
    pub kp: [i32; 8],
    pub kp_sum: i32,
    pub ap: [i32; 8],
    pub ap_avg: i32,
    /// Planetary daily character figure
    pub cp: f64,
    /// Scale cp to \[0, 9\]
    pub c9: i32,
    /// International Sunspot Number
    pub isn: i32,
    pub f10p7_obs: f64,
    pub f10p7_adj: f64,
    pub f10p7_obs_c81: f64,
    pub f10p7_obs_l81: f64,
    pub f10p7_adj_c81: f64,
    pub f10p7_adj_l81: f64,
}

fn str2num<T: core::str::FromStr>(s: &str, sidx: usize, eidx: usize) -> Option<T> {
    match s
        .chars()
        .skip(sidx)
        .take(eidx - sidx)
        .collect::<String>()
        .trim()
        .parse()
    {
        Ok(v) => Some(v),
        Err(_) => None,
    }
}

impl PartialEq for SpaceWeatherRecord {
    fn eq(&self, other: &Self) -> bool {
        self.date == other.date
    }
}

impl PartialOrd for SpaceWeatherRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.date.partial_cmp(&other.date)
    }
}

impl PartialEq<AstroTime> for SpaceWeatherRecord {
    fn eq(&self, other: &AstroTime) -> bool {
        self.date == *other
    }
}

impl PartialOrd<AstroTime> for SpaceWeatherRecord {
    fn partial_cmp(&self, other: &AstroTime) -> Option<Ordering> {
        self.date.partial_cmp(other)
    }
}

fn load_space_weather() -> SKResult<Vec<SpaceWeatherRecord>> {
    let path = datadir()
        .unwrap_or(PathBuf::from("."))
        .join("sw19571001.txt");
    if !path.is_file() {
        return skerror!("cannot load space weather data: file \"sw19571001.txt\" is missing");
    }

    let file = File::open(&path)?;

    let mut sw = Vec::<SpaceWeatherRecord>::new();

    for line in io::BufReader::new(file).lines() {
        let vline = line.unwrap();
        if vline.len() < 130 {
            continue;
        }
        if vline.chars().nth(0).unwrap() == '#' {
            continue;
        }
        let year: u32 = match str2num(&vline, 0, 4) {
            Some(v) => v,
            None => continue,
        };
        let mon: u32 = match str2num(&vline, 5, 7) {
            Some(v) => v,
            None => continue,
        };
        let day: u32 = match str2num(&vline, 8, 10) {
            Some(v) => v,
            None => continue,
        };
        sw.push(SpaceWeatherRecord {
            date: (AstroTime::from_date(year, mon, day)),
            bsrn: (str2num(&vline, 11, 15).unwrap_or(-1)),
            nd: (str2num(&vline, 16, 18).unwrap_or(-1)),
            kp: {
                let mut kparr: [i32; 8] = [-1, -1, -1, -1, -1, -1, -1, -1];
                for idx in 0..8 {
                    kparr[idx] = str2num(&vline, idx * 3 + 19, idx * 3 + 21).unwrap_or(-1);
                }
                kparr
            },
            kp_sum: (str2num(&vline, 43, 46).unwrap_or(-1)),
            ap: {
                let mut aparr: [i32; 8] = [-1, -1, -1, -1, -1, -1, -1, -1];
                for idx in 0..8 {
                    aparr[idx] = str2num(&vline, idx * 4 + 47, idx * 4 + 50).unwrap_or(-1)
                }
                aparr
            },
            ap_avg: (str2num(&vline, 79, 82).unwrap_or(-1)),
            cp: (str2num(&vline, 83, 86).unwrap_or(-1.0)),
            c9: (str2num(&vline, 87, 88).unwrap_or(-1)),
            isn: (str2num(&vline, 89, 92).unwrap_or(-1)),
            f10p7_obs: (str2num(&vline, 113, 118).unwrap_or(-1.0)),
            f10p7_adj: (str2num(&vline, 93, 98).unwrap_or(-1.0)),
            f10p7_obs_c81: (str2num(&vline, 119, 124).unwrap_or(-1.0)),
            f10p7_obs_l81: (str2num(&vline, 125, 130).unwrap_or(-1.0)),
            f10p7_adj_c81: (str2num(&vline, 101, 106).unwrap_or(-1.0)),
            f10p7_adj_l81: (str2num(&vline, 107, 112).unwrap_or(-1.0)),
        });
    } // end of for loop
    Ok(sw)
}

fn space_weather_singleton() -> &'static RwLock<SKResult<Vec<SpaceWeatherRecord>>> {
    static INSTANCE: OnceCell<RwLock<SKResult<Vec<SpaceWeatherRecord>>>> = OnceCell::new();
    INSTANCE.get_or_init(|| RwLock::new(load_space_weather()))
}

pub fn get(tm: AstroTime) -> SKResult<SpaceWeatherRecord> {
    let sw_lock = space_weather_singleton().read().unwrap();
    let sw = sw_lock.as_ref().unwrap();

    // First, try simple indexing
    let idx = (tm - sw[0].date).days().floor() as usize;
    if idx < sw.len() {
        if (tm - sw[idx].date).days().abs() < 1.0 {
            return Ok(sw[idx].clone());
        }
    }

    // More-complex lookup (is it in the future?)
    // Increase efficiency by looking backward
    let rec = sw.iter().rev().find(|x| x.date <= tm);
    if rec.is_none() {
        skerror!("Invalid date")
    } else {
        Ok(rec.unwrap().clone())
    }
}

pub fn reload() -> SKResult<()> {
    // Find writeabld data directory
    let d: Vec<PathBuf> = testdirs()
        .into_iter()
        .filter(|x| x.is_dir())
        .filter(|x| x.metadata().unwrap().permissions().readonly() == false)
        .collect();

    if d.len() == 0 {
        return skerror!("Cannot find writable data directory");
    }

    // Download most-recent EOP
    let url = "https://celestrak.org/SpaceData/sw19571001.txt";
    download_file(url, &d[0], true)?;

    *space_weather_singleton().write().unwrap() = load_space_weather();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load() {
        let tm: AstroTime = AstroTime::from_datetime(2023, 11, 14, 0, 0, 0.0);
        let r = get(tm);
        println!("r = {:?}", r);
        println!("rdate = {}", r.unwrap().date);
    }
}
