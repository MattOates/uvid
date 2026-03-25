//! Data directory resolution for reference genome files.
//!
//! Locates reference genome files on disk using platform-appropriate
//! data directories.  The search order is:
//!
//! 1. Explicit path (if the caller provides one)
//! 2. `UVID_DATA_DIR` environment variable
//! 3. Platform default via `sysdirs::data_dir() / "uvid"`
//!    - Linux:   `~/.local/share/uvid`
//!    - macOS:   `~/Library/Application Support/uvid`
//!    - Windows: `C:\Users\<user>\AppData\Roaming\uvid`
//!
//! Within the data directory, reference files are expected to be named
//! after the assembly: `GRCh38.2bit`, `GRCh37.fa`, etc.
//!
//! This module is intentionally free of UVID-specific types.

use std::path::{Path, PathBuf};

use super::error::ReferenceError;
use super::reference::{open_reference, ReferenceGenome};

/// Environment variable that overrides the default data directory.
const ENV_DATA_DIR: &str = "UVID_DATA_DIR";

/// Return the data directory for UVID reference files.
///
/// Resolution order:
/// 1. `UVID_DATA_DIR` environment variable
/// 2. `sysdirs::data_dir() / "uvid"`
///
/// Returns `None` only if no platform data directory can be determined
/// and the environment variable is not set.
pub fn data_dir() -> Option<PathBuf> {
    if let Ok(dir) = std::env::var(ENV_DATA_DIR) {
        return Some(PathBuf::from(dir));
    }
    sysdirs::data_dir().map(|d| d.join("uvid"))
}

/// File extensions to search for, in preference order.
const EXTENSIONS: &[&str] = &["2bit", "fa", "fasta"];

/// Search the data directory for a reference genome matching `assembly`.
///
/// Tries each extension in [`EXTENSIONS`] order, returning the first
/// existing file.  For example, for assembly `"GRCh38"` it checks:
///
/// - `<data_dir>/GRCh38.2bit`
/// - `<data_dir>/GRCh38.fa`
/// - `<data_dir>/GRCh38.fasta`
///
/// Returns the path to the first match, or a [`ReferenceError::NotFound`]
/// listing all paths that were tried.
pub fn find_reference(assembly: &str) -> Result<PathBuf, ReferenceError> {
    let dir = data_dir().ok_or_else(|| ReferenceError::NotFound {
        assembly: assembly.to_string(),
        searched: vec![],
    })?;

    find_reference_in(assembly, &dir)
}

/// Search a specific directory for a reference genome matching `assembly`.
pub fn find_reference_in(assembly: &str, dir: &Path) -> Result<PathBuf, ReferenceError> {
    let mut searched = Vec::new();

    for ext in EXTENSIONS {
        let candidate = dir.join(format!("{}.{}", assembly, ext));
        if candidate.exists() {
            // For FASTA, also check that the index exists
            if *ext == "fa" || *ext == "fasta" {
                let fai = {
                    let mut p = candidate.as_os_str().to_owned();
                    p.push(".fai");
                    PathBuf::from(p)
                };
                if !fai.exists() {
                    // FASTA without index — skip but record it
                    searched.push(candidate);
                    continue;
                }
            }
            return Ok(candidate);
        }
        searched.push(candidate);
    }

    Err(ReferenceError::NotFound {
        assembly: assembly.to_string(),
        searched,
    })
}

/// Open the reference genome for the given assembly name.
///
/// Searches the default data directory (or `UVID_DATA_DIR`) for a
/// matching file and opens it.
pub fn open_reference_for_assembly(
    assembly: &str,
) -> Result<Box<dyn ReferenceGenome>, ReferenceError> {
    let path = find_reference(assembly)?;
    open_reference(&path)
}

/// Open the reference genome for the given assembly, searching in a
/// specific directory.
pub fn open_reference_for_assembly_in(
    assembly: &str,
    dir: &Path,
) -> Result<Box<dyn ReferenceGenome>, ReferenceError> {
    let path = find_reference_in(assembly, dir)?;
    open_reference(&path)
}
