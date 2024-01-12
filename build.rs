use std::process::Command;

fn main() {
    let nrlmsise21_files = [
        "msis_constants.F90",
        "msis_utils.F90",
        "msis_init.F90",
        "msis_gfn.F90",
        "msis_tfn.F90",
        "msis_dfn.F90",
        "msis_calc.F90",
        "msis_gtd8d.F90",
    ];
    cc::Build::new()
        .files(
            nrlmsise21_files
                .iter()
                .map(|&file| format!("extern/nrlmsise2.1/{}", file)),
        )
        .compiler("gfortran")
        .flag("-std=legacy")
        .flag("-w")
        .flag("-O3")
        .flag("-Jextern/nrlmsise2.1/")
        .compile("nrlmsise21");

    cc::Build::new()
        .file("extern/nrlmsise/nrlmsise-00.c")
        .file("extern/nrlmsise/nrlmsise-00_data.c")
        .compile("nrlmsise");

    // Record git hash to compile-time environment variable
    let output = Command::new("git")
        .args(&["rev-parse", "HEAD"])
        .output()
        .unwrap();
    let git_hash = String::from_utf8(output.stdout).unwrap();
    println!("cargo:rustc-env=GIT_HASH={}", git_hash);
    let build_date = chrono::Utc::now().to_rfc3339();
    println!("cargo:rustc-env=BUILD_DATE={}", build_date);
    #[cfg(feature = "pybindings")]
    pyo3_build_config::add_extension_module_link_args();
}
