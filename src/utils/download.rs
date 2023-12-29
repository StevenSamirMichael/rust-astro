use crate::AstroResult;
use std::path::PathBuf;

pub fn download_file(
    url: &str,
    downloaddir: &PathBuf,
    overwrite_if_exists: bool,
) -> AstroResult<bool> {
    let fname = std::path::Path::new(url).file_name().unwrap();
    let fullpath = downloaddir.join(fname);
    if fullpath.exists() && !overwrite_if_exists {
        Ok(false)
    } else {
        println!("Downloading {}", fname.to_str().unwrap());
        let resp = ureq::get(url).call()?;

        let mut dest = std::fs::File::create(fullpath)?;
        std::io::copy(resp.into_reader().as_mut(), &mut dest)?;
        Ok(true)
    }
}
