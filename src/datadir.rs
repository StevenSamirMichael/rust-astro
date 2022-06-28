use nix::libc;
use std::fs::metadata;
use std::path::Path;
use std::{ffi::CStr, os::raw::c_void, path::PathBuf};

#[inline]
fn get_dylib_path() -> Option<PathBuf> {
    let mut dl_info = libc::Dl_info {
        dli_fname: core::ptr::null(),
        dli_fbase: core::ptr::null_mut(),
        dli_sname: core::ptr::null(),
        dli_saddr: core::ptr::null_mut(),
    };
    if unsafe {
        libc::dladdr(
            get_dylib_path as *const c_void,
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

fn find_datadir() -> Option<PathBuf> {
    fn testdir(dir: &str) {
        println!("dir = {}", dir);
    }

    let mut testdirs: Vec<PathBuf> = Vec::new();

    // Look for paths under home directory
    match std::env::var(&"HOME") {
        Ok(val) => {
            let vstr = &String::from(val);
            testdirs.push(Path::new(vstr).to_path_buf());
            testdirs.push(Path::new(vstr).join("astro-data"));
            #[cfg(target_os = "macos")]
            testdirs.push(
                Path::new(vstr)
                    .join("Libaray")
                    .join("Application Support")
                    .join("astro-data"),
            );
        }
        Err(_e) => (),
    }

    // Look for paths in environment variable
    match std::env::var(&"ASTROLIB_DATA") {
        Ok(val) => testdirs.push(Path::new(&val).to_path_buf()),
        Err(_) => (),
    }

    // Look for paths in current library directory
    match get_dylib_path() {
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

    println!("{:?}", testdirs);
    let dirname = "/Users/st16626";
    match metadata(dirname) {
        Ok(v) => {
            if v.is_dir() {
                testdir(dirname)
            }
        }
        Err(_e) => (),
    }
    Some(testdirs[0].clone())
}

static DATADIR: Option<PathBuf> = find_datadir();

pub fn get() -> Option<PathBuf> {
    DATADIR.clone()
}
