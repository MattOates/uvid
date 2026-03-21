/// 128-bit Universal Variant ID (UVID) encoding and decoding.
///
/// Bit layout (MSB to LSB):
///
/// ```text
/// Bits 127-123 (5 bits)   Chromosome       chr1-22=0-21, X=22, Y=23, M=24
/// Bits 122-91  (32 bits)  Position         1-based position within chromosome
/// Bits 90-85   (6 bits)   REF length       0-63 bases
/// Bits 84-79   (6 bits)   ALT length       0-63 bases
/// Bit  78      (1 bit)    Overflow flag    1 = sequences not stored (ref+alt > 38 bases)
/// Bits 77-2    (76 bits)  Sequence payload dynamic split [REF 2bit][ALT 2bit], zero-padded
/// Bits 1-0     (2 bits)   Assembly         0=GRCh37, 1=GRCh38, 2-3=future
/// ```
///
/// Sorting: Variants sort naturally by chromosome, then position, then by allele content.
/// Assembly is at the LSB so same-locus variants from different assemblies cluster together.
use std::fmt;

use crate::assembly::{Assembly, ChrIndex};
use crate::twobit;

// Bit field positions and masks
const CHR_SHIFT: u32 = 123;
const CHR_MASK: u128 = 0b11111;

const POS_SHIFT: u32 = 91;
const POS_MASK: u128 = 0xFFFF_FFFF; // 32 bits

const REF_LEN_SHIFT: u32 = 85;
const REF_LEN_MASK: u128 = 0b111111; // 6 bits

const ALT_LEN_SHIFT: u32 = 79;
const ALT_LEN_MASK: u128 = 0b111111; // 6 bits

const OVERFLOW_SHIFT: u32 = 78;
const OVERFLOW_MASK: u128 = 0b1;

const PAYLOAD_SHIFT: u32 = 2;
const PAYLOAD_BITS: u32 = 76;
const PAYLOAD_MASK: u128 = (1u128 << PAYLOAD_BITS) - 1;

const ASSEMBLY_MASK: u128 = 0b11;

/// Maximum number of bases that can be stored in the 76-bit payload.
pub const MAX_PAYLOAD_BASES: usize = 38; // 76 bits / 2 bits per base

/// A decoded variant representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub chr: ChrIndex,
    pub pos: u32,
    pub ref_seq: Vec<u8>,
    pub alt_seq: Vec<u8>,
    pub ref_len: usize,
    pub alt_len: usize,
    pub overflow: bool,
    pub assembly: Assembly,
}

/// A 128-bit Universal Variant ID.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Uvid128(pub u128);

impl Uvid128 {
    /// Encode a variant into a 128-bit UVID.
    ///
    /// `ref_seq` and `alt_seq` are the allele sequences as ASCII bytes.
    /// Use `b"."` or `b""` for deletion/missing alleles (encoded as length 0).
    ///
    /// Returns `None` if chromosome or position is invalid, or if allele lengths exceed 63.
    pub fn encode(
        chr: ChrIndex,
        pos: u32,
        ref_seq: &[u8],
        alt_seq: &[u8],
        assembly: Assembly,
    ) -> Option<Self> {
        // Validate chromosome index
        if chr.0 > ChrIndex::MAX {
            return None;
        }

        // Validate position (must be >= 1, fits in 32 bits)
        if pos == 0 {
            return None;
        }

        // Handle deletion symbol
        let ref_bytes = if ref_seq == b"." || ref_seq.is_empty() {
            &[] as &[u8]
        } else {
            ref_seq
        };
        let alt_bytes = if alt_seq == b"." || alt_seq.is_empty() {
            &[] as &[u8]
        } else {
            alt_seq
        };

        let ref_len = ref_bytes.len();
        let alt_len = alt_bytes.len();

        // Lengths must fit in 6 bits
        if ref_len > 63 || alt_len > 63 {
            return None;
        }

        let mut value: u128 = 0;

        // Pack chromosome (bits 127-123)
        value |= (chr.0 as u128 & CHR_MASK) << CHR_SHIFT;

        // Pack position (bits 122-91)
        value |= (pos as u128 & POS_MASK) << POS_SHIFT;

        // Pack REF length (bits 90-85)
        value |= (ref_len as u128 & REF_LEN_MASK) << REF_LEN_SHIFT;

        // Pack ALT length (bits 84-79)
        value |= (alt_len as u128 & ALT_LEN_MASK) << ALT_LEN_SHIFT;

        // Determine if sequences fit in the payload
        let total_bases = ref_len + alt_len;
        if total_bases > MAX_PAYLOAD_BASES {
            // Overflow: set flag, leave payload zeroed
            value |= OVERFLOW_MASK << OVERFLOW_SHIFT;
        } else if total_bases > 0 {
            // Pack sequences into payload
            // REF occupies the most significant portion, ALT follows
            let mut payload: u128 = 0;

            if ref_len > 0 {
                if !twobit::is_valid_dna(ref_bytes) {
                    return None;
                }
                let (ref_packed, _) = twobit::encode_sequence(ref_bytes, ref_len)?;
                // Place REF at the top of the payload
                payload |= ref_packed << (alt_len * 2);
            }

            if alt_len > 0 {
                if !twobit::is_valid_dna(alt_bytes) {
                    return None;
                }
                let (alt_packed, _) = twobit::encode_sequence(alt_bytes, alt_len)?;
                // ALT sits directly below REF in the payload
                payload |= alt_packed;
            }

            // Shift payload into position (bits 77-2)
            value |= (payload & PAYLOAD_MASK) << PAYLOAD_SHIFT;
        }

        // Pack assembly (bits 1-0)
        value |= assembly.to_bits() as u128 & ASSEMBLY_MASK;

        Some(Uvid128(value))
    }

    /// Decode a 128-bit UVID back into its component fields.
    pub fn decode(self) -> Option<Variant> {
        let value = self.0;

        let chr_bits = ((value >> CHR_SHIFT) & CHR_MASK) as u8;
        let chr = ChrIndex::from_bits(chr_bits)?;

        let pos = ((value >> POS_SHIFT) & POS_MASK) as u32;

        let ref_len = ((value >> REF_LEN_SHIFT) & REF_LEN_MASK) as usize;
        let alt_len = ((value >> ALT_LEN_SHIFT) & ALT_LEN_MASK) as usize;

        let overflow = ((value >> OVERFLOW_SHIFT) & OVERFLOW_MASK) == 1;

        let assembly_bits = (value & ASSEMBLY_MASK) as u8;
        let assembly = Assembly::from_bits(assembly_bits)?;

        let (ref_seq, alt_seq) = if overflow {
            // Sequences not stored, return empty vecs
            (Vec::new(), Vec::new())
        } else {
            let payload = (value >> PAYLOAD_SHIFT) & PAYLOAD_MASK;
            let ref_packed = if ref_len > 0 {
                payload >> (alt_len * 2)
            } else {
                0
            };
            let alt_packed = if alt_len > 0 {
                payload & ((1u128 << (alt_len * 2)) - 1)
            } else {
                0
            };
            (
                twobit::decode_sequence(ref_packed, ref_len),
                twobit::decode_sequence(alt_packed, alt_len),
            )
        };

        Some(Variant {
            chr,
            pos,
            ref_seq,
            alt_seq,
            ref_len,
            alt_len,
            overflow,
            assembly,
        })
    }

    /// Get the raw u128 value.
    pub fn as_u128(self) -> u128 {
        self.0
    }

    /// Construct from a raw u128 value.
    pub fn from_u128(value: u128) -> Self {
        Uvid128(value)
    }

    /// Compute the UVID range (lower, upper) for a genomic region query.
    ///
    /// This returns bounds such that any variant on `chr` between `start_pos`
    /// and `end_pos` (inclusive) will have a UVID in [lower, upper].
    ///
    /// All remaining fields (ref_len, alt_len, payload, assembly) are set to
    /// minimum for the lower bound and maximum for the upper bound.
    pub fn range(chr: ChrIndex, start_pos: u32, end_pos: u32) -> (Uvid128, Uvid128) {
        let chr_val = chr.0 as u128;

        // Lower bound: chr + start_pos, everything else 0
        let lower = (chr_val & CHR_MASK) << CHR_SHIFT | (start_pos as u128 & POS_MASK) << POS_SHIFT;

        // Upper bound: chr + end_pos, everything else max
        let upper = (chr_val & CHR_MASK) << CHR_SHIFT
            | (end_pos as u128 & POS_MASK) << POS_SHIFT
            | REF_LEN_MASK << REF_LEN_SHIFT
            | ALT_LEN_MASK << ALT_LEN_SHIFT
            | OVERFLOW_MASK << OVERFLOW_SHIFT
            | PAYLOAD_MASK << PAYLOAD_SHIFT
            | ASSEMBLY_MASK;

        (Uvid128(lower), Uvid128(upper))
    }

    /// Format as a hex string with dashes for readability.
    ///
    /// Format: `XXXXXXXX-XXXXXXXX-XXXXXXXX-XXXXXXXX` (32 hex digits, 4 groups of 8)
    pub fn to_hex_string(self) -> String {
        let hex = format!("{:032x}", self.0);
        format!(
            "{}-{}-{}-{}",
            &hex[0..8],
            &hex[8..16],
            &hex[16..24],
            &hex[24..32]
        )
    }

    /// Parse a UVID from a hex string (with or without dashes).
    pub fn from_hex_string(s: &str) -> Option<Self> {
        let clean: String = s.chars().filter(|c| c.is_ascii_hexdigit()).collect();
        if clean.len() > 32 {
            return None;
        }
        let value = u128::from_str_radix(&clean, 16).ok()?;
        Some(Uvid128(value))
    }
}

impl fmt::Debug for Uvid128 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Uvid128({})", self.to_hex_string())
    }
}

impl fmt::Display for Uvid128 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_hex_string())
    }
}

impl From<u128> for Uvid128 {
    fn from(value: u128) -> Self {
        Uvid128(value)
    }
}

impl From<Uvid128> for u128 {
    fn from(uvid: Uvid128) -> u128 {
        uvid.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_snp() {
        let uvid = Uvid128::encode(
            ChrIndex(2), // chr3
            1,
            b"A",
            b"G",
            Assembly::GRCh38,
        )
        .unwrap();

        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.chr, ChrIndex(2));
        assert_eq!(decoded.pos, 1);
        assert_eq!(decoded.ref_seq, b"A");
        assert_eq!(decoded.alt_seq, b"G");
        assert_eq!(decoded.ref_len, 1);
        assert_eq!(decoded.alt_len, 1);
        assert!(!decoded.overflow);
        assert_eq!(decoded.assembly, Assembly::GRCh38);
    }

    #[test]
    fn test_encode_deletion() {
        let uvid = Uvid128::encode(
            ChrIndex(0), // chr1
            100,
            b"ACGT",
            b"A",
            Assembly::GRCh38,
        )
        .unwrap();

        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.chr, ChrIndex(0));
        assert_eq!(decoded.pos, 100);
        assert_eq!(decoded.ref_seq, b"ACGT");
        assert_eq!(decoded.alt_seq, b"A");
        assert_eq!(decoded.ref_len, 4);
        assert_eq!(decoded.alt_len, 1);
    }

    #[test]
    fn test_encode_insertion() {
        let uvid = Uvid128::encode(ChrIndex(0), 500, b"A", b"ACGTACGT", Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.ref_seq, b"A");
        assert_eq!(decoded.alt_seq, b"ACGTACGT");
    }

    #[test]
    fn test_encode_dot_deletion() {
        let uvid = Uvid128::encode(ChrIndex(0), 1, b"A", b".", Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.ref_seq, b"A");
        assert_eq!(decoded.alt_seq, b"");
        assert_eq!(decoded.alt_len, 0);
    }

    #[test]
    fn test_encode_mnv() {
        let uvid = Uvid128::encode(ChrIndex(0), 1000, b"ACG", b"TGA", Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.ref_seq, b"ACG");
        assert_eq!(decoded.alt_seq, b"TGA");
    }

    #[test]
    fn test_overflow_long_sequences() {
        // 20 + 20 = 40 bases > 38 max payload
        let ref_seq = b"ACGTACGTACGTACGTACGT"; // 20 bases
        let alt_seq = b"TGCATGCATGCATGCATGCA"; // 20 bases

        let uvid = Uvid128::encode(ChrIndex(0), 1, ref_seq, alt_seq, Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert!(decoded.overflow);
        assert_eq!(decoded.ref_len, 20);
        assert_eq!(decoded.alt_len, 20);
        // Sequences are not stored in overflow mode
        assert!(decoded.ref_seq.is_empty());
        assert!(decoded.alt_seq.is_empty());
    }

    #[test]
    fn test_max_payload_boundary() {
        // Exactly 38 bases total: 19 + 19
        let ref_seq = b"ACGTACGTACGTACGTACG"; // 19 bases
        let alt_seq = b"TGCATGCATGCATGCATGC"; // 19 bases

        let uvid = Uvid128::encode(ChrIndex(0), 1, ref_seq, alt_seq, Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert!(!decoded.overflow);
        assert_eq!(decoded.ref_seq, ref_seq.to_vec());
        assert_eq!(decoded.alt_seq, alt_seq.to_vec());
    }

    #[test]
    fn test_assembly_sorting() {
        // Same variant, different assemblies -- GRCh37 should sort before GRCh38
        let uvid37 = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh37).unwrap();
        let uvid38 = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh38).unwrap();

        // Assembly is in LSB, so GRCh37 (0) sorts before GRCh38 (1)
        assert!(uvid37.0 < uvid38.0);
        // But they're close together (differ only in last 2 bits)
        assert_eq!(uvid37.0 >> 2, uvid38.0 >> 2);
    }

    #[test]
    fn test_position_sorting() {
        let uvid1 = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh38).unwrap();
        let uvid2 = Uvid128::encode(ChrIndex(0), 200, b"A", b"G", Assembly::GRCh38).unwrap();
        assert!(uvid1.0 < uvid2.0);
    }

    #[test]
    fn test_chromosome_sorting() {
        let uvid_chr1 = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh38).unwrap();
        let uvid_chr2 = Uvid128::encode(ChrIndex(1), 100, b"A", b"G", Assembly::GRCh38).unwrap();
        assert!(uvid_chr1.0 < uvid_chr2.0);
    }

    #[test]
    fn test_range_query() {
        let uvid = Uvid128::encode(ChrIndex(2), 1500, b"A", b"G", Assembly::GRCh38).unwrap();
        let (lower, upper) = Uvid128::range(ChrIndex(2), 1000, 2000);

        assert!(uvid.0 >= lower.0);
        assert!(uvid.0 <= upper.0);

        // Out of range
        let uvid_outside = Uvid128::encode(ChrIndex(2), 500, b"A", b"G", Assembly::GRCh38).unwrap();
        assert!(uvid_outside.0 < lower.0);
    }

    #[test]
    fn test_hex_string_roundtrip() {
        let uvid = Uvid128::encode(ChrIndex(2), 12345, b"GCTA", b"A", Assembly::GRCh38).unwrap();
        let hex = uvid.to_hex_string();
        let recovered = Uvid128::from_hex_string(&hex).unwrap();
        assert_eq!(uvid, recovered);
    }

    #[test]
    fn test_hex_string_no_dashes() {
        let uvid = Uvid128::encode(ChrIndex(0), 1, b"A", b"G", Assembly::GRCh38).unwrap();
        let hex = uvid.to_hex_string();
        let no_dashes = hex.replace('-', "");
        let recovered = Uvid128::from_hex_string(&no_dashes).unwrap();
        assert_eq!(uvid, recovered);
    }

    #[test]
    fn test_invalid_inputs() {
        // Position 0 is invalid
        assert!(Uvid128::encode(ChrIndex(0), 0, b"A", b"G", Assembly::GRCh38).is_none());

        // Invalid chromosome index
        assert!(Uvid128::encode(ChrIndex(25), 1, b"A", b"G", Assembly::GRCh38).is_none());

        // Invalid nucleotide
        assert!(Uvid128::encode(ChrIndex(0), 1, b"N", b"G", Assembly::GRCh38).is_none());

        // Allele too long (>63 bases)
        let long_seq = vec![b'A'; 64];
        assert!(Uvid128::encode(ChrIndex(0), 1, &long_seq, b"G", Assembly::GRCh38).is_none());
    }

    #[test]
    fn test_deterministic() {
        // Same input always produces the same UVID
        let a = Uvid128::encode(ChrIndex(2), 1, b"GCTA", b"A", Assembly::GRCh38).unwrap();
        let b = Uvid128::encode(ChrIndex(2), 1, b"GCTA", b"A", Assembly::GRCh38).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn test_different_variants_different_uvids() {
        let a = Uvid128::encode(ChrIndex(0), 1, b"A", b"G", Assembly::GRCh38).unwrap();
        let b = Uvid128::encode(ChrIndex(0), 1, b"A", b"T", Assembly::GRCh38).unwrap();
        assert_ne!(a, b);
    }
}
