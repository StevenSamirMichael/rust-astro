use std::fmt;

#[derive(Debug)]
pub struct AstroErr {
    details: String,
}

impl AstroErr {
    pub fn new(msg: &str) -> AstroErr {
        AstroErr {
            details: msg.to_string(),
        }
    }
}

impl fmt::Display for AstroErr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "AstroErr: {}", self.details)
    }
}

impl std::error::Error for AstroErr {
    fn description(&self) -> &str {
        &self.details
    }
}
