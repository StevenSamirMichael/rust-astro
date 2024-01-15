use crate::skerror;
use crate::SKResult;
use once_cell::sync::OnceCell;
use std::path::Path;
use std::{ffi::CStr, os::raw::c_void, path::PathBuf};

#[inline]
fn dylib_path() -> Option<PathBuf> {
    let mut dl_info = libc::Dl_info {
        dli_fname: core::ptr::null(),
        dli_fbase: core::ptr::null_mut(),
        dli_sname: core::ptr::null(),
        dli_saddr: core::ptr::null_mut(),
    };
    if unsafe {
        libc::dladdr(
            dylib_path as *const c_void,
            &mut dl_info as *mut libc::Dl_info,
        ) != 0
    } {
        if dl_info.dli_fname.is_null() {
            None
        } else {
            match unsafe { CStr::from_ptr(dl_info.dli_fname) }.to_str() {
                Ok(path) => Some(PathBuf::from(path)),
                Err(_) => None,
            }
        }
    } else {
        None
    }
}

pub fn testdirs() -> Vec<PathBuf> {
    let mut testdirs: Vec<PathBuf> = Vec::new();

    // Look for paths in environment variable
    match std::env::var(&"ASTROLIB_DATA") {
        Ok(val) => testdirs.push(Path::new(&val).to_path_buf()),
        Err(_) => (),
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
            testdirs.push(Path::new(vstr).join("/astro-data"));
            testdirs.push(Path::new(vstr).to_path_buf());
        }
        Err(_e) => (),
    }

    // Look for paths in current library directory
    match dylib_path() {
        Some(v) => {
            testdirs.push(Path::new(&v).join("share").join("astro-data"));
            testdirs.push(PathBuf::from(&v).join("astro-data"));
            testdirs.push(v);
        }
        None => (),
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
/// *  "ASTROLIB_DATA" environment variable
/// *  ${HOME}/share/.astro-data
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
    #[test]
    fn datadir() {
        use crate::utils::datadir;
        let d = datadir::datadir();
        assert_eq!(d.is_err(), false);
    }
}
