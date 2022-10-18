use super::settings::PropSettings;
use crate::astrotime::AstroTime;

pub trait ForceTerm<T> {
    fn ydot(&self, time: &AstroTime, state: &T) -> T;
    fn init(start_time: &AstroTime, stop_time: &AstroTime, settings: &PropSettings) -> Self
    where
        Self: Sized;
}
