pub use self::satrec::SatRec;

#[derive(PartialEq, PartialOrd, Clone, Debug, Eq, Copy)]
pub enum GravConst {
    WGS72,
    WGS72OLD,
    WGS84,
}

#[derive(PartialEq, PartialOrd, Clone, Debug, Eq, Copy)]
pub enum OpsMode {
    AFSPC,
    IMPROVED,
}

mod dpper;
mod dscom;
mod dsinit;
mod dspace;
mod getgravconst;
mod initl;
pub mod satrec;
pub mod sgp4;
mod sgp4_lowlevel;
mod sgp4init;
