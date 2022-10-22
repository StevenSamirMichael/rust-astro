pub(crate) use astroerr::astroerr;
pub use astroerr::AstroErr;
pub use astroerr::AstroResult;

mod astroerr;
pub mod datadir;
pub mod dev;

mod update_data;
pub use update_data::update_datafiles;

mod download;
pub use download::download_file;
