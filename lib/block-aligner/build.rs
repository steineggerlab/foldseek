fn main() {
    #[cfg(feature = "c_bindings")]
    {
        let crate_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
        let file = format!("{}/c/block_aligner.h", crate_dir);
        cbindgen::generate(&crate_dir)
            .unwrap()
            .write_to_file(&file);
        println!("cargo:rerun-if-changed={}", file);
    }
}
