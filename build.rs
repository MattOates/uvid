fn main() {
    // DuckDB's bundled C++ build (via duckdb-rs) uses Windows Restart Manager
    // and BCrypt APIs, but the duckdb-rs build.rs doesn't emit the necessary
    // cargo:rustc-link-lib directives for these system libraries.
    // See: https://github.com/duckdb/duckdb-rs/issues/544
    // See: https://github.com/duckdb/duckdb/blob/main/src/CMakeLists.txt
    if std::env::var("CARGO_CFG_WINDOWS").is_ok() {
        println!("cargo:rustc-link-lib=Rstrtmgr");
        println!("cargo:rustc-link-lib=Bcrypt");
    }
}
