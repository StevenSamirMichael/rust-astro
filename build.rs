fn main() {
    cc::Build::new()
        .file("extern/nrlmsise/nrlmsise-00.c")
        .file("extern/nrlmsise/nrlmsise-00_data.c")
        .compile("nrlmsise");
}
