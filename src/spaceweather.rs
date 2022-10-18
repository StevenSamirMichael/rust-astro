use std::fs::File;

use std::io::{self, BufRead};
use std::path::PathBuf;

use crate::astrotime::AstroTime;
use crate::utils::{astroerr, AstroResult, *};
use lazy_static;

#[derive(Debug)]
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
    /// Scale cp to [0,9]
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

#[derive(Debug)]
struct SpaceWeatherData {
    data: Vec<SpaceWeatherRecord>,
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

fn load_space_weather() -> AstroResult<SpaceWeatherData> {
    let path = datadir::get()
        .unwrap_or(PathBuf::from("."))
        .join("sw19571001.txt");
    if !path.is_file() {
        return astroerr!("cannot load space weather data: file \"sw19571001.txt\" is missing");
    }

    let file = File::open(&path)?;

    let mut sw = SpaceWeatherData {
        data: Vec::<SpaceWeatherRecord>::new(),
    };

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
        sw.data.push(SpaceWeatherRecord {
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

lazy_static::lazy_static! {

static ref SPACE_WEATHER_DATA: AstroResult<SpaceWeatherData> = load_space_weather();

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load() {
        let x = load_space_weather();
        println!("x = {:?}", x.unwrap().data[0]);
    }
}
