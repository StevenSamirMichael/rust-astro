pub(crate) use astroerr::astroerr;
pub use astroerr::AstroErr;
pub use astroerr::AstroResult;

mod astroerr;
mod datadir;
pub use datadir::datadir;
pub use datadir::testdirs;

#[cfg(test)]
pub mod test;

mod update_data;
pub use update_data::update_datafiles;

mod download;
pub use download::download_file;

///
/// Return git hash of compiled library
///
pub fn githash<'a>() -> &'a str {
    env!("GIT_HASH")
}

///
/// Return libary compile date
///
pub fn build_date<'a>() -> &'a str {
    env!("BUILD_DATE")
}
