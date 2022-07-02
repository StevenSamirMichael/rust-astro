//! Is this where to-level documentation goes?
//! apparently not
//! # asdfasdf

/// AstroTime
pub mod astrotime;
pub mod datadir;
pub mod earth_orientation_params;
pub mod itrfcoord;
pub mod tle;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
