use std::error::Error;
use std::fmt;

pub type AstroResult<T> = Result<T, Box<dyn Error + Send + Sync>>;

#[derive(Debug)]
pub struct AstroErr {
    details: String,
}

#[macro_export]
macro_rules! astroerr {
    ($($args:tt),*) => {{
        Err(AstroErr::new(format!($($args),*).as_str()).into())
    }};
}

pub(crate) use astroerr;

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
