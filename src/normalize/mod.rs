//! VCF variant normalization (Tan et al. 2015).
//!
//! This module implements left-alignment and parsimonious trimming of
//! VCF variant alleles so that differently-described variants that are
//! biologically identical map to the same representation.
//!
//! # Design
//!
//! The module is intentionally free of UVID-specific types (`ChrIndex`,
//! `Assembly`, etc.) and works entirely with generic primitives
//! (`&str`, `u32`, `&[u8]`).  This makes it trivially extractable into
//! a standalone crate.
//!
//! # Usage
//!
//! ```ignore
//! use uvid::normalize::{normalize, InlineReference};
//!
//! let mut reference = InlineReference::new(Default::default())
//!     .with_contig("chr1", b"ACGTACGTACGT");
//!
//! let result = normalize("chr1", 5, b"TACG", b"T", &mut reference).unwrap();
//! assert!(result.was_modified);
//! ```

pub mod algorithm;
pub mod data;
pub mod error;
pub mod reference;
#[cfg(test)]
mod tests_gatk;
#[cfg(test)]
mod tests_vt;
pub mod types;

// Re-export the most commonly used items at module level.
pub use algorithm::normalize;
pub use data::{data_dir, find_reference, open_reference_for_assembly};
pub use error::{NormalizeError, ReferenceError};
pub use reference::{
    open_reference, FastaReference, InlineReference, ReferenceGenome, TwoBitReference,
};
pub use types::NormalizedVariant;
