use nalgebra::coordinates::X;

use crate::SKResult;
use std::path::PathBuf;

pub fn download_file(
    url: &str,
    downloaddir: &PathBuf,
    overwrite_if_exists: bool,
) -> SKResult<bool> {
    let fname = std::path::Path::new(url).file_name().unwrap();
    let fullpath = downloaddir.join(fname);
    if fullpath.exists() && !overwrite_if_exists {
        Ok(false)
    } else {
        println!("Downloading {}", fname.to_str().unwrap());

        // Try to set proxy, if any, from environment variables
        let agent = ureq::AgentBuilder::new().try_proxy_from_env(true).build();

        let resp = agent.get(url).call()?;

        let mut dest = std::fs::File::create(fullpath)?;
        std::io::copy(resp.into_reader().as_mut(), &mut dest)?;
        Ok(true)
    }
}

pub fn download_to_string(url: &str) -> SKResult<String> {
    let agent = ureq::AgentBuilder::new().try_proxy_from_env(true).build();
    let resp = agent.get(url).call()?;
    let thestring = std::io::read_to_string(resp.into_reader().as_mut())?;
    Ok(thestring)
}
