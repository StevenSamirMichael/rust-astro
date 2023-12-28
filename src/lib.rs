//! Is this where to-level documentation goes?
//! # Introduction

// Type definitions
pub mod types;

/// Time and time bases (UTC, TAI, GPS, TT, etc...)
pub mod astrotime;
/// Earth orientation parameters (polar motion, delta-UT1, lenth of day)
pub mod earth_orientation_params;
/// Conversion between coordinate frames
pub mod frametransform;
/// Zonal gravity model
pub mod gravity;
/// Internation Terrestrial Reference Frame coordinates &
/// transformations to Geodetic, East-North-Up, North-East-Down
pub mod itrfcoord;
/// Solar system body ephemerides, as published by JPL
pub mod jplephem;
pub mod lpephem;
/// NRL-MISE00 Density model
pub mod nrlmsise;
/// Orbit Propagation
pub mod orbitprop;
/// SGP-4 Orbit Propagator
pub mod sgp4;
/// Solar system bodies
mod solarsystem;
pub mod spaceweather;
/// Two-line Element Set
pub mod tle;
/// Universal constants
pub mod univ;
/// Utility functions
pub mod utils;

// Integrate ordinary differential equations
pub mod ode;

// Objects available at crate level
pub use astrotime::AstroTime;
pub use astrotime::Scale as TimeScale;
pub use itrfcoord::ITRFCoord;
pub use solarsystem::SolarSystem;
pub use tle::TLE;
pub use utils::AstroErr;
pub use utils::AstroResult;

#[cfg(feature = "pybindings")]
pub mod pybindings;
