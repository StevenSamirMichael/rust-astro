extern crate reqwest;

#[cfg(test)]
type ErrResult = Box<dyn std::error::Error + Send + Sync>;

#[cfg(test)]
use super::datadir;

#[cfg(test)]
use std::io::{self, BufRead};

#[cfg(test)]
use std::path::Path;

#[cfg(test)]
use std::path::PathBuf;

#[cfg(test)]
pub fn download_if_not_exist(url: &str) -> Result<PathBuf, ErrResult> {
    let fname = match Path::new(url).file_name() {
        Some(v) => v,
        None => &std::ffi::OsStr::new("download.dat"),
    };

    let fullpath = &datadir::get().unwrap().join(fname);
    if fullpath.exists() {
        return Ok(fullpath.clone());
    }
    let resp = reqwest::blocking::get(url)?;

    let mut dest = std::fs::File::create(fullpath)?;
    let content = resp.text()?;
    std::io::copy(&mut content.as_bytes(), &mut dest)?;
    Ok(fullpath.clone())
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
#[cfg(test)]
pub fn lines_from_url(url: &str) -> Result<io::Lines<io::BufReader<std::fs::File>>, ErrResult> {
    let fullpath = download_if_not_exist(url)?;

    let file = std::fs::File::open(fullpath)?;
    let b = io::BufReader::new(file);
    Ok(b.lines())
}

#[cfg(test)]
pub fn get_project_root() -> std::io::Result<PathBuf> {
    let path = std::env::current_dir()?;
    let mut path_ancestors = path.as_path().ancestors();

    while let Some(p) = path_ancestors.next() {
        let has_cargo = std::fs::read_dir(p)?
            .into_iter()
            .any(|p| p.unwrap().file_name() == std::ffi::OsString::from("Cargo.lock"));
        if has_cargo {
            return Ok(PathBuf::from(p));
        }
    }
    Err(std::io::Error::new(
        std::io::ErrorKind::NotFound,
        "Ran out of places to find Cargo.toml",
    ))
}
