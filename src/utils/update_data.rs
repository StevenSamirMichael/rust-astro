use super::datadir::get_testdirs;
use super::download_file;
use crate::astroerr;
use crate::AstroResult;
use std::path::PathBuf;

pub fn update_datafiles(dir: Option<PathBuf>, overwrite_if_exists: bool) -> AstroResult<()> {
    // Find directory where files will be downloaded
    let downloaddir = match dir {
        Some(pb) => pb,
        None => {
            let d: Vec<PathBuf> = get_testdirs()
                .into_iter()
                .filter(|x| x.is_dir())
                .filter(|x| x.metadata().unwrap().permissions().readonly() == false)
                .collect();
            if d.len() == 0 {
                return astroerr!("Cannot find writable data directory");
            }
            d[0].clone()
        }
    };

    // List of files to download
    let urls: Vec<&str> = vec![
        "http://icgem.gfz-potsdam.de/getmodel/gfc/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.gfc",
        "http://icgem.gfz-potsdam.de/getmodel/gfc/a3375e01a717ac162962138a5e94f10466b71aa4a130d7f7d5b18ab3d5f90c3d/JGM3.gfc",
        "http://icgem.gfz-potsdam.de/getmodel/gfc/291f7d127f49fe3bcb4a633a06e8b50f461d4b13ff8e7f2046644a8148b57fc6/JGM2.gfc",
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0cbbaba92c08482d2d221c5d375264f4c6233e8f695cbcfb86ecfc5e19345a42/ITU_GRACE16.gfc",
        "https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2a.txt",
        "https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2b.txt",
        "https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2d.txt",
        "https://www.ietf.org/timezones/data/leap-seconds.list",
        "https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440",
    ];

    // List of files that are updated daily and may need to get overwritten
    let urls_overwrite: Vec<&str> = vec![
        "https://datacenter.iers.org/data/9/finals2000A.all",
        "https://celestrak.org/SpaceData/sw19571001.txt",
    ];

    // Walk through & download files
    for url in urls {
        download_file(url, &downloaddir, overwrite_if_exists)?;
    }
    // Walk through & download files
    for url in urls_overwrite {
        download_file(url, &downloaddir, true)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn update_data() {
        match update_datafiles(None, false) {
            Ok(()) => (),
            Err(e) => {
                println!("Error: {}", e.to_string());
                assert!(1 == 0);
            }
        }
    }
}
