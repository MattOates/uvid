/// 128-bit Universal Variant ID (UVID) encoding and decoding.
///
/// Bit layout (MSB to LSB):
///
/// ```text
/// Bits 127-96 (32 bits)  Linearized position  chr+pos in single coordinate space
/// Bits 95-94  (2 bits)   Assembly             0=GRCh37, 1=GRCh38, 2-3=future
/// Bit  93     (1 bit)    REF mode             0=string, 1=length
/// Bits 92-47  (46 bits)  REF payload          string: 5b len + 40b DNA; length: 46b int
/// Bit  46     (1 bit)    ALT mode             0=string, 1=length
/// Bits 45-0   (46 bits)  ALT payload          string: 5b len + 40b DNA; length: 46b int
/// ```
///
/// Sorting: Variants sort naturally by linearized position (which concatenates
/// chromosomes end-to-end), then assembly, then REF mode (string < length),
/// then REF content, then ALT mode, then ALT content.
///
/// String mode (mode bit = 0) sorts before length mode (mode bit = 1) at the
/// same locus, so variants with exact sequence encoding are grouped first.
///
/// Alleles too long (> 20 bases) or containing non-ACGT characters fall back to
/// length mode, which stores only the integer length.  Decode returns N-repeats
/// for length-mode alleles (the actual sequence can be recovered from the
/// reference genome given the position and length).
use std::fmt;

use uuid::Uuid;

use crate::allele_pack;
use crate::assembly::{self, Assembly, ChrIndex};

/// The UVID namespace UUID for UUIDv5 generation.
///
/// Computed as `UUIDv5(NAMESPACE_OID, "UVID")` where NAMESPACE_OID is the
/// ISO OID namespace from RFC 4122 (`6ba7b812-9dad-11d1-80b4-00c04fd430c8`).
///
/// Value: `2696985c-755c-53de-b6b9-1745af20d0fd`
pub const UVID_NAMESPACE: Uuid = Uuid::from_bytes([
    0x26, 0x96, 0x98, 0x5c, 0x75, 0x5c, 0x53, 0xde, 0xb6, 0xb9, 0x17, 0x45, 0xaf, 0x20, 0xd0, 0xfd,
]);

// ── Bit field positions and masks ──────────────────────────────────────────

/// Linearized position: bits 127-96 (32 bits).
const POS_SHIFT: u32 = 96;
const POS_MASK: u128 = 0xFFFF_FFFF;

/// Assembly: bits 95-94 (2 bits).
const ASM_SHIFT: u32 = 94;
const ASM_MASK: u128 = 0b11;

/// REF mode bit: bit 93.
const REF_MODE_SHIFT: u32 = 93;
const REF_MODE_MASK: u128 = 1;

/// REF payload: bits 92-47 (46 bits).
const REF_PAYLOAD_SHIFT: u32 = 47;
const REF_PAYLOAD_MASK: u128 = (1u128 << 46) - 1;

/// ALT mode bit: bit 46.
const ALT_MODE_SHIFT: u32 = 46;
const ALT_MODE_MASK: u128 = 1;

/// ALT payload: bits 45-0 (46 bits).
const ALT_PAYLOAD_MASK: u128 = (1u128 << 46) - 1;

/// A decoded variant representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub chr: ChrIndex,
    pub pos: u32,
    pub ref_seq: Vec<u8>,
    pub alt_seq: Vec<u8>,
    pub ref_len: usize,
    pub alt_len: usize,
    pub ref_is_exact: bool,
    pub alt_is_exact: bool,
    pub ref_fingerprint: Option<u32>,
    pub alt_fingerprint: Option<u32>,
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
    /// Returns `None` only if the chromosome or position is invalid (i.e. cannot
    /// be converted to a linearized coordinate).  All allele sequences are
    /// encodable: short ACGT sequences use exact 2-bit string mode, while long
    /// or non-ACGT sequences fall back to integer length mode.
    pub fn encode(
        chr: ChrIndex,
        pos: u32,
        ref_seq: &[u8],
        alt_seq: &[u8],
        assembly: Assembly,
    ) -> Option<Self> {
        // Convert chr + pos to linearized coordinate
        let linear = assembly::chr_pos_to_linear(chr, pos, assembly)?;

        // Encode each allele independently (always succeeds)
        let ref_allele = allele_pack::encode_allele(ref_seq);
        let alt_allele = allele_pack::encode_allele(alt_seq);

        // Extract mode bits (bit 45 of each 46-bit allele value)
        let ref_mode = (ref_allele >> allele_pack::MODE_BIT_POS) & 1;
        let alt_mode = (alt_allele >> allele_pack::MODE_BIT_POS) & 1;

        // Extract payloads (bits 44-0 of each 46-bit allele value)
        let ref_payload = ref_allele & allele_pack::PAYLOAD_MASK;
        let alt_payload = alt_allele & allele_pack::PAYLOAD_MASK;

        let mut value: u128 = 0;

        // Pack linearized position (bits 127-96)
        value |= (linear as u128 & POS_MASK) << POS_SHIFT;

        // Pack assembly (bits 95-94)
        value |= (assembly.to_bits() as u128 & ASM_MASK) << ASM_SHIFT;

        // Pack REF mode (bit 93)
        value |= (ref_mode as u128 & REF_MODE_MASK) << REF_MODE_SHIFT;

        // Pack REF payload (bits 92-47)
        value |= (ref_payload as u128 & REF_PAYLOAD_MASK) << REF_PAYLOAD_SHIFT;

        // Pack ALT mode (bit 46)
        value |= (alt_mode as u128 & ALT_MODE_MASK) << ALT_MODE_SHIFT;

        // Pack ALT payload (bits 45-0)
        value |= alt_payload as u128 & ALT_PAYLOAD_MASK;

        Some(Uvid128(value))
    }

    /// Decode a 128-bit UVID back into its component fields.
    ///
    /// Returns `None` if the linearized position cannot be mapped back to a
    /// valid chromosome + position, or if the assembly bits are invalid.
    pub fn decode(self) -> Option<Variant> {
        let value = self.0;

        // Extract linearized position
        let linear = ((value >> POS_SHIFT) & POS_MASK) as u32;

        // Extract assembly
        let asm_bits = ((value >> ASM_SHIFT) & ASM_MASK) as u8;
        let assembly = Assembly::from_bits(asm_bits)?;

        // Convert linearized position back to chr + pos
        let (chr, pos) = assembly::linear_to_chr_pos(linear, assembly)?;

        // Reconstruct the 46-bit allele values (mode bit + payload)
        let ref_mode = ((value >> REF_MODE_SHIFT) & REF_MODE_MASK) as u64;
        let ref_payload = ((value >> REF_PAYLOAD_SHIFT) & REF_PAYLOAD_MASK) as u64;
        let ref_allele = (ref_mode << allele_pack::MODE_BIT_POS) | ref_payload;

        let alt_mode = ((value >> ALT_MODE_SHIFT) & ALT_MODE_MASK) as u64;
        let alt_payload = (value & ALT_PAYLOAD_MASK) as u64;
        let alt_allele = (alt_mode << allele_pack::MODE_BIT_POS) | alt_payload;

        // Decode each allele
        let ref_decoded = allele_pack::decode_allele(ref_allele);
        let alt_decoded = allele_pack::decode_allele(alt_allele);

        let (ref_seq, ref_len, ref_is_exact, ref_fingerprint) = match ref_decoded {
            allele_pack::AlleleDecoded::Sequence(seq) => {
                let len = seq.len();
                (seq, len, true, None)
            }
            allele_pack::AlleleDecoded::Length { len, fingerprint } => {
                // Return N-repeats for length-mode alleles
                (vec![b'N'; len], len, false, Some(fingerprint))
            }
        };

        let (alt_seq, alt_len, alt_is_exact, alt_fingerprint) = match alt_decoded {
            allele_pack::AlleleDecoded::Sequence(seq) => {
                let len = seq.len();
                (seq, len, true, None)
            }
            allele_pack::AlleleDecoded::Length { len, fingerprint } => {
                (vec![b'N'; len], len, false, Some(fingerprint))
            }
        };

        Some(Variant {
            chr,
            pos,
            ref_seq,
            alt_seq,
            ref_len,
            alt_len,
            ref_is_exact,
            alt_is_exact,
            ref_fingerprint,
            alt_fingerprint,
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
    /// Returns bounds such that any variant on `chr` between `start_pos` and
    /// `end_pos` (inclusive) for the given assembly will have a UVID in
    /// `[lower, upper]`.
    ///
    /// Returns `None` if the chromosome/position cannot be linearized.
    pub fn range(
        chr: ChrIndex,
        start_pos: u32,
        end_pos: u32,
        assembly: Assembly,
    ) -> Option<(Uvid128, Uvid128)> {
        let linear_start = assembly::chr_pos_to_linear(chr, start_pos, assembly)?;
        let linear_end = assembly::chr_pos_to_linear(chr, end_pos, assembly)?;

        // Lower bound: position + assembly, all allele fields zero
        let lower = (linear_start as u128 & POS_MASK) << POS_SHIFT
            | (assembly.to_bits() as u128 & ASM_MASK) << ASM_SHIFT;

        // Upper bound: position + assembly, all allele fields max
        let upper = (linear_end as u128 & POS_MASK) << POS_SHIFT
            | (assembly.to_bits() as u128 & ASM_MASK) << ASM_SHIFT
            | REF_MODE_MASK << REF_MODE_SHIFT
            | REF_PAYLOAD_MASK << REF_PAYLOAD_SHIFT
            | ALT_MODE_MASK << ALT_MODE_SHIFT
            | ALT_PAYLOAD_MASK;

        Some((Uvid128(lower), Uvid128(upper)))
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

    /// Convert this UVID to a deterministic UUIDv5.
    ///
    /// Uses [`UVID_NAMESPACE`] as the namespace and the raw big-endian bytes
    /// of the 128-bit integer as the name, per RFC 4122 §4.3.
    pub fn to_uuid5(self) -> Uuid {
        Uuid::new_v5(&UVID_NAMESPACE, &self.0.to_be_bytes())
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
        assert!(decoded.ref_is_exact);
        assert!(decoded.alt_is_exact);
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
    fn test_long_sequences_use_length_mode() {
        // 21 bases each — exceeds string mode limit, falls back to length mode
        let ref_seq = b"ACGTACGTACGTACGTACGTA"; // 21 bases
        let alt_seq = b"TGCATGCATGCATGCATGCAT"; // 21 bases

        let uvid = Uvid128::encode(ChrIndex(0), 1, ref_seq, alt_seq, Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert!(!decoded.ref_is_exact);
        assert!(!decoded.alt_is_exact);
        assert_eq!(decoded.ref_len, 21);
        assert_eq!(decoded.alt_len, 21);
        // Length mode returns N-repeats
        assert_eq!(decoded.ref_seq, vec![b'N'; 21]);
        assert_eq!(decoded.alt_seq, vec![b'N'; 21]);
    }

    #[test]
    fn test_very_long_sequences_encodable() {
        // Very long sequences always succeed via length mode
        let ref_seq = vec![b'A'; 10000];
        let alt_seq = vec![b'T'; 5000];

        let uvid = Uvid128::encode(ChrIndex(0), 1, &ref_seq, &alt_seq, Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert!(!decoded.ref_is_exact);
        assert!(!decoded.alt_is_exact);
        assert_eq!(decoded.ref_len, 10000);
        assert_eq!(decoded.alt_len, 5000);
    }

    #[test]
    fn test_n_containing_sequences_use_length_mode() {
        // N bases force length mode even for short sequences
        let uvid = Uvid128::encode(ChrIndex(0), 1, b"ACGN", b"G", Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert!(!decoded.ref_is_exact); // N forced length mode
        assert!(decoded.alt_is_exact); // G is string mode
        assert_eq!(decoded.ref_len, 4);
        assert_eq!(decoded.alt_len, 1);
        assert_eq!(decoded.ref_seq, vec![b'N'; 4]);
        assert_eq!(decoded.alt_seq, b"G");
    }

    #[test]
    fn test_max_string_mode_boundary() {
        // Exactly 20 bases: string mode (max)
        let ref_seq = b"ACGTACGTACGTACGTACGT"; // 20 bases
        let alt_seq = b"TGCATGCATGCATGCATGCA"; // 20 bases

        let uvid = Uvid128::encode(ChrIndex(0), 1, ref_seq, alt_seq, Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert!(decoded.ref_is_exact);
        assert!(decoded.alt_is_exact);
        assert_eq!(decoded.ref_seq, ref_seq.to_vec());
        assert_eq!(decoded.alt_seq, alt_seq.to_vec());
    }

    #[test]
    fn test_position_sorting() {
        let uvid1 = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh38).unwrap();
        let uvid2 = Uvid128::encode(ChrIndex(0), 200, b"A", b"G", Assembly::GRCh38).unwrap();
        assert!(uvid1.0 < uvid2.0);
    }

    #[test]
    fn test_chromosome_sorting() {
        // Linearized positions: chr1 < chr2 because chr2 offset > any chr1 position
        let uvid_chr1 = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh38).unwrap();
        let uvid_chr2 = Uvid128::encode(ChrIndex(1), 100, b"A", b"G", Assembly::GRCh38).unwrap();
        assert!(uvid_chr1.0 < uvid_chr2.0);
    }

    #[test]
    fn test_assembly_sorting() {
        // Same chr+pos, different assemblies — groups by position first, then assembly
        let uvid37 = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh37).unwrap();
        let uvid38 = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh38).unwrap();

        // The linearized positions will differ between assemblies (different offset
        // tables), so sorting depends on which assembly produces a smaller linear
        // coordinate.  For chr1 pos 100 they should be the same (offset 0 for both),
        // so assembly bits at 95-94 determine order: GRCh37=0 < GRCh38=1.
        assert!(uvid37.0 < uvid38.0);
    }

    #[test]
    fn test_string_mode_sorts_before_length_mode() {
        // Same position, one allele in string mode, one in length mode
        let uvid_string = Uvid128::encode(ChrIndex(0), 1, b"ACGT", b"G", Assembly::GRCh38).unwrap();
        // Force length mode by using N-containing ref
        let uvid_length = Uvid128::encode(ChrIndex(0), 1, b"ACGN", b"G", Assembly::GRCh38).unwrap();

        // String mode (ref_mode = 0) should sort before length mode (ref_mode = 1)
        assert!(uvid_string.0 < uvid_length.0);
    }

    #[test]
    fn test_range_query() {
        let uvid = Uvid128::encode(ChrIndex(2), 1500, b"A", b"G", Assembly::GRCh38).unwrap();
        let (lower, upper) = Uvid128::range(ChrIndex(2), 1000, 2000, Assembly::GRCh38).unwrap();

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
        // Invalid chromosome index (>24)
        assert!(Uvid128::encode(ChrIndex(25), 1, b"A", b"G", Assembly::GRCh38).is_none());

        // Position 0 — chr_pos_to_linear adds (pos - 1) so pos=0 would underflow;
        // the assembly module should reject it
        assert!(Uvid128::encode(ChrIndex(0), 0, b"A", b"G", Assembly::GRCh38).is_none());
    }

    #[test]
    fn test_n_bases_no_longer_fail() {
        // N bases are handled via length mode
        let uvid = Uvid128::encode(ChrIndex(0), 1, b"N", b"G", Assembly::GRCh38);
        assert!(uvid.is_some());
    }

    #[test]
    fn test_long_allele_no_longer_fails() {
        // Alleles of any length succeed via length mode
        let long_seq = vec![b'A'; 10000];
        let uvid = Uvid128::encode(ChrIndex(0), 1, &long_seq, b"G", Assembly::GRCh38);
        assert!(uvid.is_some());
    }

    #[test]
    fn test_deterministic() {
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

    #[test]
    fn test_uuid5_deterministic() {
        let uvid = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh38).unwrap();
        let a = uvid.to_uuid5();
        let b = uvid.to_uuid5();
        assert_eq!(a, b);
    }

    #[test]
    fn test_uuid5_different_variants() {
        let a = Uvid128::encode(ChrIndex(0), 100, b"A", b"G", Assembly::GRCh38).unwrap();
        let b = Uvid128::encode(ChrIndex(0), 100, b"A", b"T", Assembly::GRCh38).unwrap();
        assert_ne!(a.to_uuid5(), b.to_uuid5());
    }

    #[test]
    fn test_uuid5_version_and_variant() {
        let uvid = Uvid128::encode(ChrIndex(2), 12345, b"GCTA", b"A", Assembly::GRCh38).unwrap();
        let uuid = uvid.to_uuid5();
        assert_eq!(uuid.get_version_num(), 5);
        assert_eq!(uuid.get_variant(), uuid::Variant::RFC4122);
    }

    #[test]
    fn test_uvid_namespace_matches_oid_derivation() {
        let expected = uuid::Uuid::new_v5(&uuid::Uuid::NAMESPACE_OID, b"UVID");
        assert_eq!(UVID_NAMESPACE, expected);
    }

    #[test]
    fn test_bit127_for_linearized_position() {
        // Linearized positions occupy the full 32 bits (genome ~3.1 billion bases).
        // Positions in later chromosomes (>~2.1B linear offset) will set bit 127,
        // which means UVIDs for these loci are negative if stored as signed i128.
        // Early chromosomes (chr1 at offset 0) have bit 127 = 0.
        let uvid_chr1 = Uvid128::encode(ChrIndex(0), 1, b"T", b"A", Assembly::GRCh38).unwrap();
        assert_eq!(uvid_chr1.0 >> 127, 0, "chr1 should have bit 127 = 0");

        // chrM end-of-genome has linear position > 2^31, so bit 127 = 1
        let uvid_m = Uvid128::encode(ChrIndex(24), 1, b"A", b"G", Assembly::GRCh38).unwrap();
        assert_eq!(
            uvid_m.0 >> 127,
            1,
            "chrM should have bit 127 = 1 (linear pos > 2^31)"
        );

        // Sorting still works correctly with u128 (unsigned): chr1 < chrM
        assert!(uvid_chr1.0 < uvid_m.0);
    }

    #[test]
    fn test_mixed_mode_alleles() {
        // REF in string mode, ALT in length mode (21 bases)
        let alt_seq = vec![b'A'; 21];
        let uvid = Uvid128::encode(ChrIndex(0), 1, b"ACG", &alt_seq, Assembly::GRCh38).unwrap();

        let decoded = uvid.decode().unwrap();
        assert!(decoded.ref_is_exact);
        assert!(!decoded.alt_is_exact);
        assert_eq!(decoded.ref_seq, b"ACG");
        assert_eq!(decoded.alt_len, 21);
    }
}
