//! Types for variant normalization results.
//!
//! This module is intentionally free of UVID-specific types.
//! All types use generic primitives so the normalize module
//! can be extracted into a standalone crate.

/// The result of normalizing a variant.
///
/// Contains the (possibly modified) chromosome, position, and allele
/// sequences after applying the
/// [Tan et al. 2015](https://doi.org/10.1093/bioinformatics/btv112)
/// normalization algorithm.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedVariant {
    /// Chromosome/contig name (unchanged from input).
    pub chrom: String,

    /// 1-based genomic position (may differ from input after left-alignment).
    pub pos: u32,

    /// Reference allele sequence, uppercase ASCII DNA bytes.
    pub ref_seq: Vec<u8>,

    /// Alternate allele sequence, uppercase ASCII DNA bytes.
    pub alt_seq: Vec<u8>,

    /// Whether the variant was modified by normalization.
    ///
    /// `false` if the variant was already in normalized form, or if
    /// normalization was skipped (e.g. SNV, symbolic allele).
    pub was_modified: bool,
}

impl NormalizedVariant {
    /// Create a new normalized variant that was not modified.
    ///
    /// Convenience constructor for the pass-through case (SNVs,
    /// symbolic alleles, already-normalized variants).
    pub fn unchanged(chrom: &str, pos: u32, ref_seq: &[u8], alt_seq: &[u8]) -> Self {
        NormalizedVariant {
            chrom: chrom.to_string(),
            pos,
            ref_seq: ref_seq.to_vec(),
            alt_seq: alt_seq.to_vec(),
            was_modified: false,
        }
    }
}
