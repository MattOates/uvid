//! Two-bit nucleotide encoding using conventional bioinformatics ordering.
//!
//! Encoding: A=0b00, C=0b01, G=0b10, T=0b11
//!
//! This gives a clean 4-element lookup for decode and a simple match for encode.
//!
//! ## Per-allele 46-bit encoding
//!
//! Each allele (REF or ALT) is encoded into 46 bits using a mode bit:
//!
//! ```text
//! Mode 0 (string mode, bit 45 = 0): exact 2-bit encoded sequence
//!   Bit 45:      0 (string mode)
//!   Bits 44-40:  5-bit sequence length (0-20)
//!   Bits 39-0:   up to 20 bases of 2-bit DNA (40 bits), zero-padded right
//!
//! Mode 1 (length mode, bit 45 = 1): length + Rabin fingerprint
//!   Bit 45:      1 (length mode)
//!   Bits 44-17:  28-bit integer length (max 268,435,455 — covers chr1 at 248M)
//!   Bits 16-0:   17-bit Rabin fingerprint of the sequence
//! ```
//!
//! String mode sorts before length mode (0 < 1) at the same locus.
//! In length mode, the actual sequence is not stored; decode returns N-repeats.
//!
//! ## Rabin fingerprint for sequence discrimination
//!
//! Length-mode alleles store a 17-bit Rabin fingerprint of the sequence to
//! reduce collisions between variants at the same locus with the same allele
//! length but different sequences (e.g. insertions of nearly-identical tandem
//! repeats differing by a few SNPs).
//!
//! The fingerprint uses polynomial division in GF(2) with an irreducible
//! polynomial of degree 17.  The DNA sequence is treated as a bitstream where
//! each base contributes 2 bits (A=00, C=01, G=10, T=11), shifted into the
//! register MSB-first.  Non-ACGT characters (e.g. N) are treated as A (0b00),
//! which is the most common base; this is acceptable because N-containing
//! sequences in length mode are extremely rare and the fingerprint is a
//! best-effort discriminator, not a guarantee.
//!
//! **Irreducible polynomial**: `x^17 + x^3 + 1` (binary: `0x20009`).
//! This is a known irreducible polynomial over GF(2) chosen for its sparsity
//! (only 3 non-zero terms), which minimises the cost of the reduction step.
//!
//! **Properties**:
//! - Position-sensitive: "ACG" and "GCA" produce different fingerprints
//!   because the polynomial shift encodes positional information.
//! - Single-base perturbation sensitive: changing one base anywhere in the
//!   sequence changes the remainder of the polynomial division, producing
//!   a different fingerprint with probability ~(1 - 2^-17) ≈ 99.999%.
//! - Deterministic and cheap: only bitwise shifts and XORs, no multiplication.
//! - 1-in-131,072 false collision rate for random same-length sequences.

/// Encode a single nucleotide character to its 2-bit representation.
///
/// Returns `None` for non-ACGT characters.
#[inline]
pub fn encode_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0b00),
        b'C' | b'c' => Some(0b01),
        b'G' | b'g' => Some(0b10),
        b'T' | b't' => Some(0b11),
        _ => None,
    }
}

/// Decode a 2-bit value back to a nucleotide ASCII byte.
///
/// Only the lowest 2 bits of `bits` are considered.
#[inline]
pub fn decode_base(bits: u8) -> u8 {
    const LOOKUP: [u8; 4] = [b'A', b'C', b'G', b'T'];
    LOOKUP[(bits & 0b11) as usize]
}

/// Encode a DNA sequence into a u128, packed from the MSB side.
///
/// Returns the packed bits and the number of bases encoded.
/// Bases are packed left-to-right starting at the most significant position
/// within the allocated bit range.
///
/// `max_bases` limits how many bases will be encoded (for fitting into a bit budget).
/// Returns `None` if any character is not a valid nucleotide.
pub fn encode_sequence(seq: &[u8], max_bases: usize) -> Option<(u128, usize)> {
    let count = seq.len().min(max_bases);
    let mut packed: u128 = 0;
    for (i, &base) in seq[..count].iter().enumerate() {
        let bits = encode_base(base)? as u128;
        // Pack from the left: first base goes into the highest bit positions
        // within a 76-bit payload, the caller is responsible for shifting into
        // the correct position within the full 128-bit word.
        let shift = (count - 1 - i) * 2;
        packed |= bits << shift;
    }
    Some((packed, count))
}

/// Decode a 2-bit packed sequence back into a byte vector of nucleotides.
///
/// `packed` contains the bases packed from MSB side.
/// `count` is the number of bases to decode.
pub fn decode_sequence(packed: u128, count: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(count);
    for i in 0..count {
        let shift = (count - 1 - i) * 2;
        let bits = ((packed >> shift) & 0b11) as u8;
        seq.push(decode_base(bits));
    }
    seq
}

/// Check if a byte slice is a valid DNA sequence (only A, C, G, T characters).
pub fn is_valid_dna(seq: &[u8]) -> bool {
    seq.iter()
        .all(|&b| matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
}

// ── Per-allele 46-bit mode encoding ────────────────────────────────────────

/// Number of bits in each allele field.
pub const ALLELE_BITS: u32 = 46;

/// Maximum number of bases that fit in string mode (40 payload bits / 2 bits per base).
pub const MAX_STRING_BASES: usize = 20;

/// Mask for the full 46-bit allele field.
pub const ALLELE_MASK: u64 = (1u64 << ALLELE_BITS) - 1;

/// The mode bit position within a 46-bit allele field.
///
/// Bit 45 = mode: 0 for string, 1 for length.
/// Public so that `uvid128` can split allele values into mode + payload.
pub const MODE_BIT_POS: u32 = 45;

/// Mask for the mode bit.
const MODE_MASK: u64 = 1u64 << MODE_BIT_POS;

/// The mode bit position (internal alias).
const MODE_BIT: u32 = MODE_BIT_POS;

/// Mask for bits below the mode bit (the 45-bit payload).
///
/// Public so that `uvid128` can extract the payload from an allele value.
pub const PAYLOAD_MASK: u64 = MODE_MASK - 1; // bits 44-0

/// Number of bits for the sequence length prefix in string mode.
const STRING_LEN_BITS: u32 = 5;

/// Shift for the string length field (below mode bit, above sequence payload).
const STRING_LEN_SHIFT: u32 = MODE_BIT - STRING_LEN_BITS; // 40

/// Mask for the string length field (5 bits).
const STRING_LEN_MASK: u64 = ((1u64 << STRING_LEN_BITS) - 1) << STRING_LEN_SHIFT;

// (STRING_PAYLOAD_BITS removed — value equals STRING_LEN_SHIFT and was unused)

// ── Length mode sub-layout constants ───────────────────────────────────────
//
// Length mode (bit 45 = 1): bits 44-0 split into 28-bit length + 17-bit Rabin hash.
//
// ```text
// Bit 45:      1 (length mode)
// Bits 44-17:  28-bit integer length (max 268,435,455)
// Bits 16-0:   17-bit Rabin fingerprint
// ```

/// Number of bits for the Rabin fingerprint in length mode.
pub const FINGERPRINT_BITS: u32 = 17;

/// Number of bits for the integer length in length mode.
pub const LENGTH_BITS: u32 = 28;

/// Shift for the length field within the 45-bit payload (above the fingerprint).
const LENGTH_SHIFT: u32 = FINGERPRINT_BITS; // 17

/// Mask for the 17-bit fingerprint (bits 16-0 of the payload).
const FINGERPRINT_MASK: u64 = (1u64 << FINGERPRINT_BITS) - 1;

/// Mask for the 28-bit length field (bits 44-17 of the payload).
const LENGTH_FIELD_MASK: u64 = ((1u64 << LENGTH_BITS) - 1) << LENGTH_SHIFT;

/// Maximum length encodable in length mode (2^28 - 1 = 268,435,455).
pub const MAX_LENGTH_MODE_LEN: usize = (1usize << LENGTH_BITS) - 1;

/// Irreducible polynomial for Rabin fingerprint: x^17 + x^3 + 1 (binary).
///
/// Only the lower 17 bits are used for XOR reduction (the leading x^17 term
/// is implicit in the shift-and-reduce algorithm).
const RABIN_POLY: u32 = 0x0_0009; // bits: x^3 + x^0 = 0b0000_0000_0000_1001

/// Result of decoding a 46-bit allele field.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlleleDecoded {
    /// String mode: the exact DNA sequence was stored.
    Sequence(Vec<u8>),
    /// Length mode: the integer length and Rabin fingerprint were stored.
    /// The actual bases are unknown (decode as N-repeats).
    Length { len: usize, fingerprint: u32 },
}

/// Compute a 17-bit Rabin fingerprint of a DNA sequence.
///
/// The sequence is treated as a bitstream where each base contributes 2 bits
/// (A=00, C=01, G=10, T=11), shifted into a polynomial division register
/// MSB-first. Non-ACGT characters (including N) are treated as A (0b00).
///
/// Uses irreducible polynomial x^17 + x^3 + 1 over GF(2).
///
/// The result is in the range `[0, 2^17)` (0 to 131,071).
pub fn rabin_fingerprint_dna(seq: &[u8]) -> u32 {
    let mut hash: u32 = 0;

    for &base in seq {
        let bits: u32 = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => 0b00, // N and other non-ACGT → A
        };

        // Shift in 2 bits at a time (MSB first within the 2-bit pair).
        // We need to check for overflow after each single-bit shift.

        // Shift in high bit of the 2-bit pair
        hash <<= 1;
        hash |= (bits >> 1) & 1;
        if hash & (1 << FINGERPRINT_BITS) != 0 {
            hash ^= RABIN_POLY;
            hash &= FINGERPRINT_MASK as u32;
        }

        // Shift in low bit of the 2-bit pair
        hash <<= 1;
        hash |= bits & 1;
        if hash & (1 << FINGERPRINT_BITS) != 0 {
            hash ^= RABIN_POLY;
            hash &= FINGERPRINT_MASK as u32;
        }
    }

    hash
}

/// Encode a single allele (REF or ALT) into a 46-bit packed value.
///
/// - Empty or `"."` → string mode with length 0 (represents deletion/missing).
/// - Length 1-20 with valid ACGT bases → string mode with exact 2-bit sequence.
/// - Length > 20 OR contains non-ACGT bases → length mode with integer length.
///
/// This function always succeeds (no `None` return).
pub fn encode_allele(seq: &[u8]) -> u64 {
    // Handle deletion symbol or empty
    if seq.is_empty() || seq == b"." {
        // String mode, length 0, no payload
        return 0;
    }

    let len = seq.len();

    // Try string mode: must be <= 20 bases and all valid ACGT
    if len <= MAX_STRING_BASES && is_valid_dna(seq) {
        // String mode: bit 45 = 0
        let mut packed: u64 = 0;

        // Pack length into bits 44-40
        packed |= (len as u64) << STRING_LEN_SHIFT;

        // Pack 2-bit encoded bases into bits 39-0 (from MSB side)
        for (i, &base) in seq.iter().enumerate() {
            let bits = encode_base(base).unwrap() as u64;
            // First base at highest position (bits 39-38), etc.
            let shift = (MAX_STRING_BASES - 1 - i) * 2;
            packed |= bits << shift;
        }

        packed
    } else {
        // Length mode: bit 45 = 1, bits 44-17 = 28-bit length, bits 16-0 = 17-bit Rabin hash
        let fingerprint = rabin_fingerprint_dna(seq) as u64;
        MODE_MASK | ((len as u64) << LENGTH_SHIFT) | (fingerprint & FINGERPRINT_MASK)
    }
}

/// Decode a 46-bit allele field back to its components.
///
/// Returns `AlleleDecoded::Sequence` for string mode (exact bases) or
/// `AlleleDecoded::Length` for length mode (integer length only).
pub fn decode_allele(packed: u64) -> AlleleDecoded {
    if packed & MODE_MASK != 0 {
        // Length mode: bits 44-17 are the 28-bit length, bits 16-0 are the 17-bit fingerprint
        let payload = packed & PAYLOAD_MASK;
        let len = ((payload & LENGTH_FIELD_MASK) >> LENGTH_SHIFT) as usize;
        let fingerprint = (payload & FINGERPRINT_MASK) as u32;
        AlleleDecoded::Length { len, fingerprint }
    } else {
        // String mode: bits 44-40 are the length, bits 39-0 are 2-bit DNA
        let len = ((packed & STRING_LEN_MASK) >> STRING_LEN_SHIFT) as usize;
        if len == 0 {
            return AlleleDecoded::Sequence(Vec::new());
        }
        let mut seq = Vec::with_capacity(len);
        for i in 0..len {
            let shift = (MAX_STRING_BASES - 1 - i) * 2;
            let bits = ((packed >> shift) & 0b11) as u8;
            seq.push(decode_base(bits));
        }
        AlleleDecoded::Sequence(seq)
    }
}

/// Check if a 46-bit allele field is in string mode (exact sequence stored).
#[inline]
pub fn is_string_mode(packed: u64) -> bool {
    packed & MODE_MASK == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A'), Some(0b00));
        assert_eq!(encode_base(b'C'), Some(0b01));
        assert_eq!(encode_base(b'G'), Some(0b10));
        assert_eq!(encode_base(b'T'), Some(0b11));
        // Case insensitive
        assert_eq!(encode_base(b'a'), Some(0b00));
        assert_eq!(encode_base(b'c'), Some(0b01));
        assert_eq!(encode_base(b'g'), Some(0b10));
        assert_eq!(encode_base(b't'), Some(0b11));
        // Invalid
        assert_eq!(encode_base(b'N'), None);
        assert_eq!(encode_base(b'X'), None);
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(0b00), b'A');
        assert_eq!(decode_base(0b01), b'C');
        assert_eq!(decode_base(0b10), b'G');
        assert_eq!(decode_base(0b11), b'T');
    }

    #[test]
    fn test_roundtrip_base() {
        for base in [b'A', b'C', b'G', b'T'] {
            let encoded = encode_base(base).unwrap();
            let decoded = decode_base(encoded);
            assert_eq!(decoded, base);
        }
    }

    #[test]
    fn test_encode_sequence() {
        // Single base
        let (packed, count) = encode_sequence(b"A", 38).unwrap();
        assert_eq!(count, 1);
        assert_eq!(packed, 0b00);

        let (packed, count) = encode_sequence(b"T", 38).unwrap();
        assert_eq!(count, 1);
        assert_eq!(packed, 0b11);

        // GCTA
        let (packed, count) = encode_sequence(b"GCTA", 38).unwrap();
        assert_eq!(count, 4);
        // G=10, C=01, T=11, A=00 -> 10_01_11_00 = 0x9C = 156
        assert_eq!(packed, 0b10_01_11_00);

        // Max bases limit
        let (packed, count) = encode_sequence(b"GCTAGCTA", 4).unwrap();
        assert_eq!(count, 4);
        assert_eq!(packed, 0b10_01_11_00);
    }

    #[test]
    fn test_decode_sequence() {
        let seq = decode_sequence(0b10_01_11_00, 4);
        assert_eq!(&seq, b"GCTA");

        let seq = decode_sequence(0b00, 1);
        assert_eq!(&seq, b"A");

        let seq = decode_sequence(0b11, 1);
        assert_eq!(&seq, b"T");
    }

    #[test]
    fn test_roundtrip_sequence() {
        let cases = [b"A".as_slice(), b"GCTA", b"ACGTACGTACGTACGT", b"TTTTTTTTTT"];
        for case in cases {
            let (packed, count) = encode_sequence(case, 38).unwrap();
            let decoded = decode_sequence(packed, count);
            assert_eq!(decoded, case);
        }
    }

    #[test]
    fn test_encode_sequence_invalid() {
        assert!(encode_sequence(b"ACGN", 38).is_none());
        assert!(encode_sequence(b"HELLO", 38).is_none());
    }

    #[test]
    fn test_is_valid_dna() {
        assert!(is_valid_dna(b"ACGT"));
        assert!(is_valid_dna(b"acgt"));
        assert!(is_valid_dna(b"AcGt"));
        assert!(!is_valid_dna(b"ACGN"));
        // Empty is vacuously valid (`.all()` returns true for empty iterators)
        assert!(is_valid_dna(b""));
    }

    #[test]
    fn test_long_sequence_truncation() {
        let long_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 44 bases
        let (packed, count) = encode_sequence(long_seq, 38).unwrap();
        assert_eq!(count, 38);
        let decoded = decode_sequence(packed, count);
        assert_eq!(&decoded, &long_seq[..38]);
    }

    // ── Per-allele 46-bit encoding tests ────────────────────────────────

    #[test]
    fn test_encode_allele_empty() {
        // Empty sequence → string mode, length 0
        let packed = encode_allele(b"");
        assert_eq!(packed, 0);
        assert!(is_string_mode(packed));
        assert_eq!(decode_allele(packed), AlleleDecoded::Sequence(Vec::new()));
    }

    #[test]
    fn test_encode_allele_dot() {
        // Deletion symbol "." → string mode, length 0
        let packed = encode_allele(b".");
        assert_eq!(packed, 0);
        assert!(is_string_mode(packed));
        assert_eq!(decode_allele(packed), AlleleDecoded::Sequence(Vec::new()));
    }

    #[test]
    fn test_encode_allele_single_base() {
        // Single base "A" → string mode, length 1
        let packed = encode_allele(b"A");
        assert!(is_string_mode(packed));
        // Length = 1 in bits 44-40, A=00 at bit position 39-38
        let expected_len = 1u64 << STRING_LEN_SHIFT;
        let expected_base = 0b00u64 << ((MAX_STRING_BASES - 1) * 2); // A=00 at MSB
        assert_eq!(packed, expected_len | expected_base);
        assert_eq!(
            decode_allele(packed),
            AlleleDecoded::Sequence(b"A".to_vec())
        );
    }

    #[test]
    fn test_encode_allele_snv() {
        // "G" → string mode
        let packed = encode_allele(b"G");
        assert!(is_string_mode(packed));
        assert_eq!(
            decode_allele(packed),
            AlleleDecoded::Sequence(b"G".to_vec())
        );
    }

    #[test]
    fn test_encode_allele_short_sequence() {
        // "GCTA" → string mode, length 4
        let packed = encode_allele(b"GCTA");
        assert!(is_string_mode(packed));
        assert_eq!(
            decode_allele(packed),
            AlleleDecoded::Sequence(b"GCTA".to_vec())
        );
    }

    #[test]
    fn test_encode_allele_max_string() {
        // 20 bases → exactly at string mode limit
        let seq = b"ACGTACGTACGTACGTACGT"; // 20 bases
        let packed = encode_allele(seq);
        assert!(is_string_mode(packed));
        assert_eq!(decode_allele(packed), AlleleDecoded::Sequence(seq.to_vec()));
    }

    #[test]
    fn test_encode_allele_overflow_to_length() {
        // 21 bases → length mode
        let seq = b"ACGTACGTACGTACGTACGTA"; // 21 bases
        let packed = encode_allele(seq);
        assert!(!is_string_mode(packed));
        let decoded = decode_allele(packed);
        match decoded {
            AlleleDecoded::Length { len, fingerprint } => {
                assert_eq!(len, 21);
                assert_eq!(fingerprint, rabin_fingerprint_dna(seq));
            }
            _ => panic!("Expected Length variant"),
        }
    }

    #[test]
    fn test_encode_allele_long_sequence() {
        // Very long sequence → length mode
        let seq = vec![b'A'; 10000];
        let packed = encode_allele(&seq);
        assert!(!is_string_mode(packed));
        match decode_allele(packed) {
            AlleleDecoded::Length { len, fingerprint } => {
                assert_eq!(len, 10000);
                assert_eq!(fingerprint, rabin_fingerprint_dna(&seq));
            }
            _ => panic!("Expected Length variant"),
        }
    }

    #[test]
    fn test_encode_allele_n_bases() {
        // N-containing sequence falls back to length mode even if short
        let packed = encode_allele(b"ACGN");
        assert!(!is_string_mode(packed));
        match decode_allele(packed) {
            AlleleDecoded::Length { len, fingerprint } => {
                assert_eq!(len, 4);
                assert_eq!(fingerprint, rabin_fingerprint_dna(b"ACGN"));
            }
            _ => panic!("Expected Length variant"),
        }
    }

    #[test]
    fn test_encode_allele_case_insensitive() {
        let upper = encode_allele(b"ACGT");
        let lower = encode_allele(b"acgt");
        assert_eq!(upper, lower);
    }

    #[test]
    fn test_allele_roundtrip_string_mode() {
        let cases: &[&[u8]] = &[
            b"",
            b"A",
            b"C",
            b"G",
            b"T",
            b"AT",
            b"GCTA",
            b"ACGTACGT",
            b"ACGTACGTACGTACGTACGT", // 20 bases (max)
        ];
        for case in cases {
            let packed = encode_allele(case);
            assert!(
                is_string_mode(packed),
                "Expected string mode for {:?}",
                std::str::from_utf8(case)
            );
            match decode_allele(packed) {
                AlleleDecoded::Sequence(seq) => assert_eq!(&seq, case),
                AlleleDecoded::Length { .. } => {
                    panic!("Expected Sequence for {:?}", std::str::from_utf8(case))
                }
            }
        }
    }

    #[test]
    fn test_allele_roundtrip_length_mode() {
        let cases = [21, 63, 100, 1000, 9983, 100000];
        for len in cases {
            let seq = vec![b'A'; len];
            let packed = encode_allele(&seq);
            assert!(
                !is_string_mode(packed),
                "Expected length mode for len={}",
                len
            );
            match decode_allele(packed) {
                AlleleDecoded::Length {
                    len: decoded_len,
                    fingerprint,
                } => {
                    assert_eq!(decoded_len, len);
                    assert_eq!(fingerprint, rabin_fingerprint_dna(&seq));
                }
                _ => panic!("Expected Length variant for len={}", len),
            }
        }
    }

    #[test]
    fn test_allele_string_sorts_before_length() {
        // String mode (mode bit 0) should sort before length mode (mode bit 1)
        let string_packed = encode_allele(b"ACGT"); // string mode
        let length_packed = encode_allele(b"ACGTACGTACGTACGTACGTACGT"); // 24 bases, length mode
        assert!(string_packed < length_packed);
    }

    #[test]
    fn test_allele_fits_46_bits() {
        // Verify no allele encoding exceeds 46 bits
        let cases: &[&[u8]] = &[
            b"",
            b"A",
            b"ACGTACGTACGTACGTACGT",  // string mode
            b"ACGTACGTACGTACGTACGTA", // length mode, 21
        ];
        for case in cases {
            let packed = encode_allele(case);
            assert_eq!(
                packed & !ALLELE_MASK,
                0,
                "Exceeds 46 bits for {:?}",
                std::str::from_utf8(case)
            );
        }
        // Max length mode value: mode=1, length=MAX_LENGTH_MODE_LEN, fingerprint=all 1s
        let packed = MODE_MASK | ((MAX_LENGTH_MODE_LEN as u64) << LENGTH_SHIFT) | FINGERPRINT_MASK;
        assert_eq!(packed & !ALLELE_MASK, 0, "Max length mode exceeds 46 bits");
    }

    #[test]
    fn test_allele_deterministic() {
        let packed1 = encode_allele(b"GCTA");
        let packed2 = encode_allele(b"GCTA");
        assert_eq!(packed1, packed2);
    }

    // ── Rabin fingerprint tests ─────────────────────────────────────────

    #[test]
    fn test_rabin_empty() {
        assert_eq!(rabin_fingerprint_dna(b""), 0);
    }

    #[test]
    fn test_rabin_single_bases() {
        // Single bases should produce small nonzero values (except A=00)
        assert_eq!(rabin_fingerprint_dna(b"A"), 0); // A=00, shifting in zeros
        assert_ne!(rabin_fingerprint_dna(b"C"), 0);
        assert_ne!(rabin_fingerprint_dna(b"G"), 0);
        assert_ne!(rabin_fingerprint_dna(b"T"), 0);
    }

    #[test]
    fn test_rabin_deterministic() {
        let a = rabin_fingerprint_dna(b"ACGTACGTACGTACGTACGTACGT");
        let b = rabin_fingerprint_dna(b"ACGTACGTACGTACGTACGTACGT");
        assert_eq!(a, b);
    }

    #[test]
    fn test_rabin_position_sensitive() {
        // Different orderings of the same bases must produce different fingerprints
        let a = rabin_fingerprint_dna(b"ACG");
        let b = rabin_fingerprint_dna(b"GCA");
        assert_ne!(a, b);
    }

    #[test]
    fn test_rabin_single_snp_differs() {
        // Sequences differing by one base should (almost certainly) differ
        let a = rabin_fingerprint_dna(b"ACGTACGTACGTACGTACGTACGT");
        let b = rabin_fingerprint_dna(b"ACGTACGTACGTACGTACGTACGG");
        assert_ne!(a, b, "Single SNP should change fingerprint");
    }

    #[test]
    fn test_rabin_fits_17_bits() {
        // Fingerprint must be < 2^17 = 131072
        let seqs: &[&[u8]] = &[
            b"A",
            b"ACGTACGTACGT",
            b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            b"GCTAGCTAGCTAGCTAGCTAGCTAGCTA",
        ];
        for seq in seqs {
            let fp = rabin_fingerprint_dna(seq);
            assert!(
                fp < (1 << FINGERPRINT_BITS),
                "Fingerprint {} exceeds {} bits for {:?}",
                fp,
                FINGERPRINT_BITS,
                std::str::from_utf8(seq)
            );
        }
    }

    #[test]
    fn test_rabin_n_treated_as_a() {
        // N is treated as A (0b00)
        let with_n = rabin_fingerprint_dna(b"NCGT");
        let with_a = rabin_fingerprint_dna(b"ACGT");
        assert_eq!(with_n, with_a);
    }

    #[test]
    fn test_rabin_case_insensitive() {
        let upper = rabin_fingerprint_dna(b"ACGTACGT");
        let lower = rabin_fingerprint_dna(b"acgtacgt");
        assert_eq!(upper, lower);
    }

    #[test]
    fn test_length_mode_different_sequences_differ() {
        // Same length, different sequences → same length in UVID, but different fingerprints
        let seq_a = b"ACGTACGTACGTACGTACGTACGT"; // 24 bases
        let seq_b = b"ACGTACGTACGTACGTACGTACGG"; // 24 bases, last base differs
        let packed_a = encode_allele(seq_a);
        let packed_b = encode_allele(seq_b);
        assert_ne!(
            packed_a, packed_b,
            "Different sequences should produce different packed values"
        );
    }

    #[test]
    fn test_length_mode_roundtrip_with_fingerprint() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGT"; // 28 bases
        let packed = encode_allele(seq);
        assert!(!is_string_mode(packed));

        let expected_fp = rabin_fingerprint_dna(seq);
        match decode_allele(packed) {
            AlleleDecoded::Length { len, fingerprint } => {
                assert_eq!(len, 28);
                assert_eq!(fingerprint, expected_fp);
            }
            _ => panic!("Expected Length variant"),
        }
    }
}
