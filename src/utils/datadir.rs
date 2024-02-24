use crate::skerror;
use crate::SKResult;
use once_cell::sync::OnceCell;
use process_path::get_dylib_path;
use std::path::Path;
use std::path::PathBuf;

pub fn testdirs() -> Vec<PathBuf> {
    let mut testdirs: Vec<PathBuf> = Vec::new();

    // Look for paths in environment variable
    match std::env::var(&"ASTRO_DATA") {
        Ok(val) => testdirs.push(Path::new(&val).to_path_buf()),
        Err(_) => (),
    }

    // Look for paths in current library directory
    match get_dylib_path() {
        Some(v) => {
            #[cfg(feature = "pybindings")]
            testdirs.push(Path::new(&v).parent().unwrap().join("astro-data"));
            testdirs.push(
                Path::new(&v)
                    .parent()
                    .unwrap()
                    .join("share")
                    .join("astro-data"),
            );
            testdirs.push(v);
        }
        None => (),
    }
    // Look for paths under home directory
    match std::env::var(&"HOME") {
        Ok(val) => {
            let vstr = &String::from(val);

            #[cfg(target_os = "macos")]
            testdirs.push(
                Path::new(vstr)
                    .join("Library")
                    .join("Application Support")
                    .join("astro-data"),
            );
            testdirs.push(Path::new(vstr).join("astro-data"));
            testdirs.push(Path::new(vstr).to_path_buf());
        }
        Err(_e) => (),
    }

    testdirs.push(Path::new(&"/usr/share/astrodata").to_path_buf());

    // On mac, look in root library directory
    #[cfg(target_os = "macos")]
    testdirs.push(Path::new(&"/Library/Application Support/astro-data").to_path_buf());

    testdirs
}

/// Get directory where astronomy data is stored
///
/// Tries the following paths in order, and stops when the
/// files are found
///
/// *  "ASTRO_DATA" environment variable
/// *  ${HOME}/astro-data
/// *  ${HOME}
/// *  /usr/share/astro-data
/// *  On Mac Only:
///    * /Library/Application Support/astro-data
///    * ${Home}/Library/Application Support/astro-data
///
/// Returns:
///
///  * Option<<std::path::PathBuf>> representing directory
///    where files are stored
///
pub fn datadir() -> SKResult<PathBuf> {
    static INSTANCE: OnceCell<SKResult<PathBuf>> = OnceCell::new();
    let res = INSTANCE.get_or_init(|| {
        for ref dir in testdirs() {
            let p = PathBuf::from(&dir).join("tab5.2a.txt");
            if p.is_file() {
                return Ok(dir.to_path_buf().clone());
            }
        }
        skerror!("Could not find valid data directory.")
    });
    match res.as_ref() {
        Ok(v) => Ok(v.clone()),
        Err(_e) => skerror!("Could not find valid data directory."),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dylib() {
        let p = get_dylib_path();
        println!("p = {:?}", p);
    }

    #[test]
    fn datadir() {
        use crate::utils::datadir;
        let d = datadir::datadir();
        assert_eq!(d.is_err(), false);
    }
}
