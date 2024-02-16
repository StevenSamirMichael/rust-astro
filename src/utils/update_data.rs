use super::download_file;
use super::testdirs;
use crate::skerror;
use crate::SKResult;
use json::JsonValue;
use std::path::PathBuf;

pub fn update_datafiles(dir: Option<PathBuf>, overwrite_if_exists: bool) -> SKResult<()> {
    // Find directory where files will be downloaded
    let downloaddir = match dir {
        Some(pb) => pb,
        None => {
            let d: Vec<PathBuf> = testdirs()
                .into_iter()
                .filter(|x| x.is_dir())
                .filter(|x| x.metadata().unwrap().permissions().readonly() == false)
                .collect();
            if d.len() == 0 {
                return skerror!("Cannot find writable data directory");
            }
            d[0].clone()
        }
    };
    println!("Downloading to {}", downloaddir.to_str().unwrap());

    // List of files to download
    let urls: Vec<&str> = vec![
        "http://icgem.gfz-potsdam.de/getmodel/gfc/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.gfc",
        "http://icgem.gfz-potsdam.de/getmodel/gfc/a3375e01a717ac162962138a5e94f10466b71aa4a130d7f7d5b18ab3d5f90c3d/JGM3.gfc",
        "http://icgem.gfz-potsdam.de/getmodel/gfc/291f7d127f49fe3bcb4a633a06e8b50f461d4b13ff8e7f2046644a8148b57fc6/JGM2.gfc",
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0cbbaba92c08482d2d221c5d375264f4c6233e8f695cbcfb86ecfc5e19345a42/ITU_GRACE16.gfc",
        "https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2a.txt",
        "https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2b.txt",
        "https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2d.txt",
        "https://ftp.iana.org/tz/tzdb-2020a/leap-seconds.list",
        "https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440",
    ];

    // List of files that are updated daily and may need to get overwritten
    let urls_overwrite: Vec<&str> = vec![
        "https://datacenter.iers.org/data/9/finals2000A.all",
        "https://celestrak.org/SpaceData/sw19571001.txt",
    ];

    let mut joiners: Vec<std::thread::JoinHandle<SKResult<bool>>> = Vec::new();

    // Walk through & download files
    for url in urls {
        let d = downloaddir.clone();
        joiners.push(std::thread::spawn(move || {
            download_file(url, &d, overwrite_if_exists.clone())
        }));
    }
    // Walk through & download files
    for url in urls_overwrite {
        let d = downloaddir.clone();
        joiners.push(std::thread::spawn(move || download_file(url, &d, true)));
    }

    // Wait for all the threads to funish
    for jh in joiners {
        jh.join().unwrap()?;
    }

    Ok(())
}

fn download_from_json(v: &JsonValue, basedir: std::path::PathBuf, baseurl: String) -> SKResult<()> {
    if v.is_object() {
        let r1: Vec<SKResult<()>> = v
            .entries()
            .map(|entry: (&str, &JsonValue)| -> SKResult<()> {
                let pbnew = basedir.join(entry.0);
                if !pbnew.is_dir() {
                    std::fs::create_dir_all(pbnew.clone())?;
                }
                let mut newurl = baseurl.clone();
                newurl.push_str(format!("/{}", entry.0).as_str());
                download_from_json(entry.1, pbnew.clone(), newurl)?;
                Ok(())
            })
            .filter(|res| match res {
                Ok(_) => false,
                Err(_) => true,
            })
            .collect();
        if r1.len() > 0 {
            return skerror!("Could not parse entries");
        }
    } else if v.is_array() {
        let r2: Vec<SKResult<()>> = v
            .members()
            .map(|val| -> SKResult<()> {
                download_from_json(val, basedir.clone(), baseurl.clone())?;
                Ok(())
            })
            .filter(|res| match res {
                Ok(_) => false,
                Err(_) => true,
            })
            .collect();
        if r2.len() > 0 {
            return skerror!("could not parse array entries");
        }
    } else if v.is_string() {
        let mut newurl = baseurl.clone();
        newurl.push_str(format!("/{}", v).as_str());
        println!("Downloading {} into {:?}", newurl, basedir);
        download_file(newurl.as_str(), &basedir, false)?;
    } else {
        return skerror!("invalid json for downloading files??!!");
    }

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn parse_json() {
        let sj = std::fs::read_to_string("satkit-testvecs/files.json").unwrap();
        let fmap = json::parse(sj.as_str()).unwrap();
        match download_from_json(
            &fmap,
            std::path::PathBuf::new().join("testdl"),
            String::from("https://stevensamirmichael.github.io/satkit-testvecs/"),
        ) {
            Ok(()) => {}
            Err(e) => {
                println!("Error: {}", e.to_string());
                assert!(1 == 0);
            }
        }
        println!("hi steven");
    }

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
