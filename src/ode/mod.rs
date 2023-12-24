mod nalgebra;
pub mod rk_adaptive;
pub mod rk_adaptive_settings;
pub mod rk_explicit;
//mod rkf45;
mod rkts54;
mod rkv65;
mod rkv65_table;
mod rkv87;
mod rkv87_table;
mod rkv98;
mod rkv98_table;
pub mod types;

pub use rk_adaptive::RKAdaptive;
pub use rk_adaptive_settings::RKAdaptiveSettings;

pub mod solvers {
    pub use super::rkts54::RKTS54;
    pub use super::rkv65::RKV65;
    pub use super::rkv87::RKV87;
    pub use super::rkv98::RKV98;
}

pub use types::*;

#[cfg(test)]
mod harmonic_oscillator;
#[cfg(test)]
use harmonic_oscillator::HarmonicOscillator;
