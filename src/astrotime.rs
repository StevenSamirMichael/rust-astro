#[derive(PartialEq, PartialOrd, Copy, Clone)]
pub struct AstroTime {
    mjd: f64,
}

impl AstroTime {
    #[inline]
    pub fn new(mjd: f64) -> AstroTime {
        AstroTime { mjd: mjd }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn datadir() {
        assert_eq!(1, 1);
    }
}
