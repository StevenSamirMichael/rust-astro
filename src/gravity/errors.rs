use std::fmt;

#[derive(Debug)]
pub struct GravityError {
    details: String,
}

impl GravityError {
    pub fn new(msg: &str) -> GravityError {
        GravityError {
            details: msg.to_string(),
        }
    }
}

impl fmt::Display for GravityError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.details)
    }
}

impl std::error::Error for GravityError {
    fn description(&self) -> &str {
        &self.details
    }
}
