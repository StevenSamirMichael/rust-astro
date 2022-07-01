pub mod astrotime;
pub mod datadir;
pub mod earth_orientation_params;
pub mod itrfcoord;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
