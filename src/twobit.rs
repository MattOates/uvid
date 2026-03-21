/// Two-bit nucleotide encoding using conventional bioinformatics ordering.
///
/// Encoding: A=0b00, C=0b01, G=0b10, T=0b11
///
/// This gives a clean 4-element lookup for decode and a simple match for encode.

/// Encode a single nucleotide character to its 2-bit representation.
///
/// Returns `None` for invalid characters.
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
}
