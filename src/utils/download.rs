use crate::AstroResult;
use reqwest;
use std::path::PathBuf;

pub fn download_file(
    url: &str,
    downloaddir: &PathBuf,
    overwrite_if_exists: bool,
) -> AstroResult<bool> {
    let fname = std::path::Path::new(url).file_name().unwrap();
    let fullpath = downloaddir.join(fname);
    if fullpath.exists() && !overwrite_if_exists {
        println!(
            "File: {} already exists.  Not downloading",
            fullpath.as_os_str().to_str().unwrap()
        );
        Ok(false)
    } else {
        println!("Downloading {}", fname.to_str().unwrap());
        let resp = reqwest::blocking::get(url)?;

        let mut dest = std::fs::File::create(fullpath)?;
        let content = resp.text()?;
        std::io::copy(&mut content.as_bytes(), &mut dest)?;
        Ok(true)
    }
}
