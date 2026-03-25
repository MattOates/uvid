//! Error types for variant normalization.
//!
//! This module is intentionally free of UVID-specific types.
//! All errors use generic primitives (strings, positions) so the
//! normalize module can be extracted into a standalone crate.

use std::fmt;
use std::io;
use std::path::PathBuf;

/// Errors that can occur when accessing a reference genome.
#[derive(Debug)]
pub enum ReferenceError {
    /// The requested chromosome/contig was not found in the reference.
    ChromosomeNotFound(String),

    /// The requested range is out of bounds for the chromosome.
    OutOfBounds {
        chrom: String,
        requested_start: u32,
        requested_end: u32,
        chrom_length: u32,
    },

    /// No reference file was found for the given assembly.
    NotFound {
        assembly: String,
        searched: Vec<PathBuf>,
    },

    /// The reference file format is not supported.
    UnsupportedFormat(PathBuf),

    /// An I/O error occurred while reading the reference.
    Io(io::Error),
}

impl fmt::Display for ReferenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReferenceError::ChromosomeNotFound(chrom) => {
                write!(f, "chromosome not found in reference: {}", chrom)
            }
            ReferenceError::OutOfBounds {
                chrom,
                requested_start,
                requested_end,
                chrom_length,
            } => {
                write!(
                    f,
                    "position {}:{}-{} out of bounds (chromosome length: {})",
                    chrom, requested_start, requested_end, chrom_length
                )
            }
            ReferenceError::NotFound { assembly, searched } => {
                write!(
                    f,
                    "no reference genome found for assembly '{}'; searched: {:?}",
                    assembly, searched
                )
            }
            ReferenceError::UnsupportedFormat(path) => {
                write!(f, "unsupported reference format: {}", path.display())
            }
            ReferenceError::Io(e) => write!(f, "reference I/O error: {}", e),
        }
    }
}

impl std::error::Error for ReferenceError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ReferenceError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for ReferenceError {
    fn from(e: io::Error) -> Self {
        ReferenceError::Io(e)
    }
}

/// Errors that can occur during variant normalization.
#[derive(Debug)]
pub enum NormalizeError {
    /// A reference genome error occurred during normalization.
    Reference(ReferenceError),

    /// Normalization did not converge within the iteration limit.
    ///
    /// This is a safety guard — in practice, convergence should always
    /// happen within a few iterations for well-formed variants.
    DidNotConverge {
        chrom: String,
        pos: u32,
        iterations: u32,
    },
}

impl fmt::Display for NormalizeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NormalizeError::Reference(e) => write!(f, "normalization failed: {}", e),
            NormalizeError::DidNotConverge {
                chrom,
                pos,
                iterations,
            } => {
                write!(
                    f,
                    "normalization did not converge for {}:{} after {} iterations",
                    chrom, pos, iterations
                )
            }
        }
    }
}

impl std::error::Error for NormalizeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            NormalizeError::Reference(e) => Some(e),
            _ => None,
        }
    }
}

impl From<ReferenceError> for NormalizeError {
    fn from(e: ReferenceError) -> Self {
        NormalizeError::Reference(e)
    }
}
