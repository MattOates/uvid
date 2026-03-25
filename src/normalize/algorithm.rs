//! Variant normalization algorithm
//! ([Tan et al. 2015](https://doi.org/10.1093/bioinformatics/btv112)).
//!
//! Implements the "unified representation" normalization:
//! 1. **Right-trim**: strip matching trailing bases from REF and ALT
//! 2. **Left-trim**: strip matching leading bases, adjusting POS
//! 3. **Left-align**: while last bases match, prepend upstream reference
//!    base and shift POS left
//!
//! Steps 1-3 repeat until convergence (no changes in an iteration).
//!
//! This module is intentionally free of UVID-specific types.

use super::error::NormalizeError;
use super::reference::ReferenceGenome;
use super::types::NormalizedVariant;

/// Maximum iterations before we declare non-convergence.
/// In practice, normalization converges in 1-2 iterations.
const MAX_ITERATIONS: u32 = 100;

/// Returns `true` if the allele is symbolic (e.g. `<DEL>`, `<INS>`)
/// or a breakend (`]` / `[` notation), or missing (`.`), or the
/// spanning deletion (`*`).
fn is_symbolic_or_special(allele: &[u8]) -> bool {
    if allele.is_empty() {
        return true;
    }
    // Missing allele
    if allele == b"." {
        return true;
    }
    // Star (spanning deletion)
    if allele == b"*" {
        return true;
    }
    // Symbolic: <...>
    if allele.first() == Some(&b'<') {
        return true;
    }
    // Breakend notation
    if allele.contains(&b'[') || allele.contains(&b']') {
        return true;
    }
    false
}

/// Normalize a single variant using the
/// [Tan et al. 2015](https://doi.org/10.1093/bioinformatics/btv112) algorithm.
///
/// # Arguments
///
/// * `chrom` — Chromosome/contig name (e.g. `"chr1"`, `"20"`)
/// * `pos` — 1-based VCF position
/// * `ref_seq` — Reference allele bytes (uppercase ASCII DNA)
/// * `alt_seq` — Alternate allele bytes (uppercase ASCII DNA)
/// * `reference` — Reference genome for fetching upstream bases
///
/// # Returns
///
/// A [`NormalizedVariant`] with the (possibly adjusted) position and
/// allele sequences.  If the variant was already normalized or was
/// skipped (SNV, symbolic), `was_modified` is `false`.
///
/// # Skipped variant types
///
/// - Symbolic alleles (`<DEL>`, `<INS>`, etc.)
/// - Missing alleles (`.`)
/// - Spanning deletions (`*`)
/// - Breakend notation
pub fn normalize(
    chrom: &str,
    pos: u32,
    ref_seq: &[u8],
    alt_seq: &[u8],
    reference: &mut dyn ReferenceGenome,
) -> Result<NormalizedVariant, NormalizeError> {
    // Skip symbolic/special alleles
    if is_symbolic_or_special(ref_seq) || is_symbolic_or_special(alt_seq) {
        return Ok(NormalizedVariant::unchanged(chrom, pos, ref_seq, alt_seq));
    }

    // Skip if either allele is empty (shouldn't happen in valid VCF, but be safe)
    if ref_seq.is_empty() || alt_seq.is_empty() {
        return Ok(NormalizedVariant::unchanged(chrom, pos, ref_seq, alt_seq));
    }

    // Work with mutable copies
    let mut r: Vec<u8> = ref_seq.to_vec();
    let mut a: Vec<u8> = alt_seq.to_vec();
    let mut p = pos;
    let orig_pos = pos;

    for _iter in 0..MAX_ITERATIONS {
        let (r2, a2, p2) = normalize_one_pass(chrom, p, &r, &a, reference)?;

        if r2 == r && a2 == a && p2 == p {
            // Converged
            let was_modified = p != orig_pos || r.as_slice() != ref_seq || a.as_slice() != alt_seq;
            return Ok(NormalizedVariant {
                chrom: chrom.to_string(),
                pos: p,
                ref_seq: r,
                alt_seq: a,
                was_modified,
            });
        }

        r = r2;
        a = a2;
        p = p2;
    }

    Err(NormalizeError::DidNotConverge {
        chrom: chrom.to_string(),
        pos: orig_pos,
        iterations: MAX_ITERATIONS,
    })
}

/// One pass of the normalization algorithm: right-trim, left-trim, left-align.
fn normalize_one_pass(
    chrom: &str,
    mut pos: u32,
    ref_seq: &[u8],
    alt_seq: &[u8],
    reference: &mut dyn ReferenceGenome,
) -> Result<(Vec<u8>, Vec<u8>, u32), NormalizeError> {
    let mut r = ref_seq.to_vec();
    let mut a = alt_seq.to_vec();

    // --- Step 1: Right-trim ---
    // Remove matching bases from the end, keeping at least 1 base each.
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }

    // --- Step 2: Left-trim ---
    // Remove matching bases from the start, adjusting POS.
    // Keep at least 1 base each.
    while r.len() > 1 && a.len() > 1 && r.first() == a.first() {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }

    // --- Step 3: Left-align ---
    // Only for indels: one allele is a prefix of the other, or after
    // trimming we have one allele of length 1.
    // While the last bases of REF and ALT match, prepend an upstream
    // reference base and shift POS left.
    //
    // Skip left-alignment for complex/MNV variants (both > 1 base
    // after trimming, and neither is a pure prefix of the other).
    let is_indel = r.len() == 1 || a.len() == 1;

    if is_indel {
        while pos > 1 && r.last() == a.last() {
            // Shift left: prepend the base at pos-2 (0-based: pos-1 is current,
            // so upstream is pos-2)
            let upstream_pos = pos - 1; // 1-based position of the base to prepend
            let upstream_base = reference.get_sequence(chrom, upstream_pos - 1, upstream_pos)?; // 0-based [pos-2, pos-1)

            if upstream_base.is_empty() {
                break;
            }

            let base = upstream_base[0];

            // Remove last base (they match)
            r.pop();
            a.pop();

            // Prepend upstream base
            r.insert(0, base);
            a.insert(0, base);

            pos -= 1;
        }
    }

    Ok((r, a, pos))
}

// ===========================================================================
// Unit tests — derived from bcftools norm test suite (MIT license)
// ===========================================================================
//
// Source: https://github.com/samtools/bcftools, test/norm.{vcf,out,fa}
// License: MIT
//
// Each test case records its origin for attribution.

#[cfg(test)]
mod tests {
    use super::super::reference::InlineReference;
    use super::*;
    use std::collections::HashMap;

    /// Test case for normalization.
    struct NormTestCase {
        name: &'static str,
        source: &'static str,
        chrom: &'static str,
        pos: u32,
        ref_seq: &'static [u8],
        alt_seq: &'static [u8],
        expected_pos: u32,
        expected_ref: &'static [u8],
        expected_alt: &'static [u8],
    }

    /// Build the inline reference matching bcftools test/norm.fa.
    fn bcftools_norm_reference() -> InlineReference {
        InlineReference::new(HashMap::new())
            // >20 20:1339000-1339300
            .with_contig(
                "20",
                b"AGGATGGGGCTCATTAATAGAGCTCCACTTGTCTCCAGAATCACTGGTGAGGAAGGGGAG\
                  TGTTGCCCCCACATTCGTGCACAGCAGGGATGGTTCACCGAACTCCACACCAGTCTCTGC\
                  AGAGCCTGTTGGGGAGAGGAGGGCTGTGGTTTCTTTGATGGTGTTCACCTGGAGTAGAGC\
                  AAGTATTGTCAAAAGGGTCATCCTCGGAGGTTGCAGTGAGCCGAGATCGCACCATTGCAC\
                  TGCAGCCTGGGAGACAGAGCAAGACTCCATCTCAAAAAAAAAAAAAAAAAAAAAGGCCAT\
                  C",
            )
            // >1 1:10143-10443
            .with_contig(
                "1",
                b"CTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACC\
                  CTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTA\
                  ACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA\
                  ACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTA\
                  ACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC",
            )
            // >2 1:1382388-1382602
            .with_contig(
                "2",
                b"GGGCGTCTCATAGCTGGAGCAATGGCGAGCGCCTGGACAAGGGAGGGGAAGGGGTTCTTA\
                  TTACTGACGCGGGTAGCCCCTACTGCTGTGTGGTTCCCCTATTTTTTTTTTTTTCTTTTT\
                  GAGACGGAGTCTCGCTCTGTCACCCAGGCTGGAGTGCAGTGGCACAATCTCGGCTCACTG\
                  CAAGCTCCACCTCCTGGGTTCACGCCATTCTCCTG",
            )
            // >3 madeup
            .with_contig(
                "3",
                b"ACTGGACACGTGGACACACACACACACACACACACACACACAGTCAAACCACCTACCAGA",
            )
            // >4 20:8917026-8917085
            .with_contig(
                "4",
                b"TCCCCTCTTGACCTCTCTCTATTTTTTTTTTTTTTTCTGAGATGGATTTTTGCTCTTGTT",
            )
            // >5 20:18724313-18724343
            .with_contig("5", b"GTCTCAAAAAAAAAAAAAAAAAAAAGAAAAG")
            // >21
            .with_contig(
                "21",
                b"TTTATTATTATTATTATTAAATTGAATTTATTTAGTGTACATACATTCATGTGTATTGTG",
            )
            // >22
            .with_contig("22", b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
    }

    fn run_test_cases(cases: &[NormTestCase]) {
        let mut reference = bcftools_norm_reference();
        for tc in cases {
            let result = normalize(tc.chrom, tc.pos, tc.ref_seq, tc.alt_seq, &mut reference)
                .unwrap_or_else(|e| panic!("{}: normalization failed: {}", tc.name, e));
            assert_eq!(
                result.pos, tc.expected_pos,
                "{} ({}): pos mismatch: got {}, expected {}",
                tc.name, tc.source, result.pos, tc.expected_pos
            );
            assert_eq!(
                result.ref_seq,
                tc.expected_ref,
                "{} ({}): ref mismatch: got {:?}, expected {:?}",
                tc.name,
                tc.source,
                String::from_utf8_lossy(&result.ref_seq),
                String::from_utf8_lossy(tc.expected_ref)
            );
            assert_eq!(
                result.alt_seq,
                tc.expected_alt,
                "{} ({}): alt mismatch: got {:?}, expected {:?}",
                tc.name,
                tc.source,
                String::from_utf8_lossy(&result.alt_seq),
                String::from_utf8_lossy(tc.expected_alt)
            );
        }
    }

    // ------------------------------------------------------------------
    // SNV: pass-through (no normalization needed)
    // ------------------------------------------------------------------

    #[test]
    fn test_snv_passthrough() {
        run_test_cases(&[
            NormTestCase {
                name: "snv_passthrough_chr20",
                source: "bcftools:norm.vcf (MIT, github.com/samtools/bcftools)",
                chrom: "20",
                pos: 81,
                ref_seq: b"A",
                alt_seq: b"C",
                expected_pos: 81,
                expected_ref: b"A",
                expected_alt: b"C",
            },
            NormTestCase {
                name: "snv_passthrough_chr2",
                source: "bcftools:norm.vcf (MIT, github.com/samtools/bcftools)",
                chrom: "2",
                pos: 115,
                ref_seq: b"C",
                alt_seq: b"T",
                expected_pos: 115,
                expected_ref: b"C",
                expected_alt: b"T",
            },
        ]);
    }

    // ------------------------------------------------------------------
    // Right-trim only
    // ------------------------------------------------------------------

    #[test]
    fn test_right_trim_overspecified_snv() {
        // 20:80 CACAG -> CACAT becomes G -> T at pos 84
        run_test_cases(&[NormTestCase {
            name: "right_trim_overspecified_snv",
            source: "bcftools:norm.vcf:20:80 (MIT, github.com/samtools/bcftools)",
            chrom: "20",
            pos: 80,
            ref_seq: b"CACAG",
            alt_seq: b"CACAT",
            expected_pos: 84,
            expected_ref: b"G",
            expected_alt: b"T",
        }]);
    }

    #[test]
    fn test_trim_to_snv() {
        // 20:95 TCACCG -> ACACCG becomes T -> A at pos 95
        run_test_cases(&[NormTestCase {
            name: "trim_to_snv",
            source: "bcftools:norm.vcf:20:95 TCACCG/ACACCG (MIT, github.com/samtools/bcftools)",
            chrom: "20",
            pos: 95,
            ref_seq: b"TCACCG",
            alt_seq: b"ACACCG",
            expected_pos: 95,
            expected_ref: b"T",
            expected_alt: b"A",
        }]);
    }

    // ------------------------------------------------------------------
    // Complex / MNV: trim only (no left-alignment)
    // ------------------------------------------------------------------

    #[test]
    fn test_complex_mnv_trim_only() {
        // 20:95 TCACCG -> AAAAAA stays as-is (complex, both > 1 base)
        run_test_cases(&[NormTestCase {
            name: "complex_mnv_no_leftalign",
            source: "bcftools:norm.vcf:20:95 TCACCG/AAAAAA (MIT, github.com/samtools/bcftools)",
            chrom: "20",
            pos: 95,
            ref_seq: b"TCACCG",
            alt_seq: b"AAAAAA",
            expected_pos: 95,
            expected_ref: b"TCACCG",
            expected_alt: b"AAAAAA",
        }]);
    }

    #[test]
    fn test_mnv_trim_ref_and_alt() {
        // 20:3 GATG -> GACT: right-trim G, left-trim G => pos 5 TG -> CT
        // (MNV after trim: both 2 bases, no left-alignment)
        run_test_cases(&[NormTestCase {
            name: "mnv_trim_gatg_gact",
            source: "bcftools:norm.vcf:20:3 GATG/GACT (MIT, github.com/samtools/bcftools)",
            chrom: "20",
            pos: 3,
            ref_seq: b"GATG",
            alt_seq: b"GACT",
            expected_pos: 5,
            expected_ref: b"TG",
            expected_alt: b"CT",
        }]);
    }

    // ------------------------------------------------------------------
    // Left-align through homopolymer
    // ------------------------------------------------------------------

    #[test]
    fn test_left_align_poly_t_insertion() {
        // 2:101 ATTTTTTTTTTTTT -> ATTTTTTTTTTTTTTT becomes A -> ATT at pos 101
        run_test_cases(&[NormTestCase {
            name: "left_align_poly_t_insertion",
            source: "bcftools:norm.vcf:2:101 (MIT, github.com/samtools/bcftools)",
            chrom: "2",
            pos: 101,
            ref_seq: b"ATTTTTTTTTTTTT",
            alt_seq: b"ATTTTTTTTTTTTTTT",
            expected_pos: 101,
            expected_ref: b"A",
            expected_alt: b"ATT",
        }]);
    }

    #[test]
    fn test_left_align_indel_chr2_114() {
        // 2:114 TC -> TTCC becomes T -> TTC at pos 114
        run_test_cases(&[NormTestCase {
            name: "left_align_indel_chr2_114_alt1",
            source: "bcftools:norm.vcf:2:114 TC/TTCC (MIT, github.com/samtools/bcftools)",
            chrom: "2",
            pos: 114,
            ref_seq: b"TC",
            alt_seq: b"TTCC",
            expected_pos: 114,
            expected_ref: b"T",
            expected_alt: b"TTC",
        }]);
    }

    #[test]
    fn test_left_align_indel_chr2_114_alt2() {
        // 2:114 TC -> TTC: right-trim C gives T -> TT at pos 114.
        // Per-allele left-alignment shifts through the poly-T run (pos 102-114)
        // to anchor on the A at pos 101, yielding A -> AT at pos 101.
        //
        // Note: bcftools normalizes multiallelic records as a unit, so it
        // keeps this at pos 114.  Our per-allele normalization produces a
        // different (equally correct) representation.
        run_test_cases(&[NormTestCase {
            name: "left_align_indel_chr2_114_alt2",
            source: "bcftools:norm.vcf:2:114 TC/TTC (MIT, github.com/samtools/bcftools)",
            chrom: "2",
            pos: 114,
            ref_seq: b"TC",
            alt_seq: b"TTC",
            expected_pos: 101,
            expected_ref: b"A",
            expected_alt: b"AT",
        }]);
    }

    // ------------------------------------------------------------------
    // Left-align in poly-A run
    // ------------------------------------------------------------------

    #[test]
    fn test_left_align_poly_a_insertion_chr20_273() {
        // 20:273 CAAAAAAAAAAAAAAAAAAAAA -> CAAAAAAAAAAAAAAAAAAAAAAA
        // becomes C -> CAA at pos 273
        run_test_cases(&[NormTestCase {
            name: "left_align_poly_a_ins_2bp",
            source: "bcftools:norm.vcf:20:273 21A/23A (MIT, github.com/samtools/bcftools)",
            chrom: "20",
            pos: 273,
            ref_seq: b"CAAAAAAAAAAAAAAAAAAAAA",
            alt_seq: b"CAAAAAAAAAAAAAAAAAAAAAAA",
            expected_pos: 273,
            expected_ref: b"C",
            expected_alt: b"CAA",
        }]);
    }

    #[test]
    fn test_left_align_poly_a_insertion_chr20_273_4bp() {
        // 20:273 CAAAAAAAAAAAAAAAAAAAAA -> CAAAAAAAAAAAAAAAAAAAAAAAA
        // becomes C -> CAAA at pos 273
        run_test_cases(&[NormTestCase {
            name: "left_align_poly_a_ins_4bp",
            source: "bcftools:norm.vcf:20:273 21A/25A (MIT, github.com/samtools/bcftools)",
            chrom: "20",
            pos: 273,
            ref_seq: b"CAAAAAAAAAAAAAAAAAAAAA",
            alt_seq: b"CAAAAAAAAAAAAAAAAAAAAAAAA",
            expected_pos: 273,
            expected_ref: b"C",
            expected_alt: b"CAAA",
        }]);
    }

    #[test]
    fn test_left_align_poly_a_insertion_chr20_274() {
        // 20:274 AAAAAAAAA -> AAAAAAAAAAAAAAAAAAA
        // becomes C -> CAAAAAAAAAA at pos 273
        run_test_cases(&[NormTestCase {
            name: "left_align_poly_a_ins_from_274",
            source: "bcftools:norm.vcf:20:274 (MIT, github.com/samtools/bcftools)",
            chrom: "20",
            pos: 274,
            ref_seq: b"AAAAAAAAA",
            alt_seq: b"AAAAAAAAAAAAAAAAAAA",
            expected_pos: 273,
            expected_ref: b"C",
            expected_alt: b"CAAAAAAAAAA",
        }]);
    }

    #[test]
    fn test_left_align_poly_a_insertion_chr20_278() {
        // 20:278 AAAAAAAAAAAAAAAAA -> AAAAAAAAAAAAAAAAAAA
        // becomes C -> CAA at pos 273
        run_test_cases(&[NormTestCase {
            name: "left_align_poly_a_ins_from_278",
            source: "bcftools:norm.vcf:20:278 (MIT, github.com/samtools/bcftools)",
            chrom: "20",
            pos: 278,
            ref_seq: b"AAAAAAAAAAAAAAAAA",
            alt_seq: b"AAAAAAAAAAAAAAAAAAA",
            expected_pos: 273,
            expected_ref: b"C",
            expected_alt: b"CAA",
        }]);
    }

    // ------------------------------------------------------------------
    // Dinucleotide repeat left-alignment
    // ------------------------------------------------------------------

    #[test]
    fn test_left_align_ca_repeat_deletion() {
        // 3:15 CACA -> CAC becomes CA -> C at pos 17
        run_test_cases(&[NormTestCase {
            name: "left_align_ca_repeat_del",
            source: "bcftools:norm.vcf:3:15 CACA/CAC (MIT, github.com/samtools/bcftools)",
            chrom: "3",
            pos: 15,
            ref_seq: b"CACA",
            alt_seq: b"CAC",
            expected_pos: 17,
            expected_ref: b"CA",
            expected_alt: b"C",
        }]);
    }

    // ------------------------------------------------------------------
    // Left-align GA insertion on chr5
    // ------------------------------------------------------------------

    #[test]
    fn test_left_align_ga_insertion() {
        // 5:22 A -> AGA becomes A -> AAG at pos 21
        run_test_cases(&[NormTestCase {
            name: "left_align_ga_insertion",
            source: "bcftools:norm.vcf:5:22 A/AGA (MIT, github.com/samtools/bcftools)",
            chrom: "5",
            pos: 22,
            ref_seq: b"A",
            alt_seq: b"AGA",
            expected_pos: 21,
            expected_ref: b"A",
            expected_alt: b"AAG",
        }]);
    }

    // ------------------------------------------------------------------
    // Large complex variant: right-trim only
    // ------------------------------------------------------------------

    #[test]
    fn test_large_complex_right_trim() {
        // 2:1 long REF -> ACGT: right-trim trailing T from both
        // REF loses last T: ...CCACC, ALT loses last T: ACG, pos stays 1
        run_test_cases(&[NormTestCase {
            name: "large_complex_right_trim",
            source: "bcftools:norm.vcf:2:1 (MIT, github.com/samtools/bcftools)",
            chrom: "2",
            pos: 1,
            ref_seq: b"GGGCGTCTCATAGCTGGAGCAATGGCGAGCGCCTGGACAAGGGAGGGGAAGGGGTTCTTATTACTGACGCGGGTAGCCCCTACTGCTGTGTGGTTCCCCTATTTTTTTTTTTTTCTTTTTGAGACGGAGTCTCGCTCTGTCACCCAGGCTGGAGTGCAGTGGCACAATCTCGGCTCACTGCAAGCTCCACCT",
            alt_seq: b"ACGT",
            expected_pos: 1,
            expected_ref: b"GGGCGTCTCATAGCTGGAGCAATGGCGAGCGCCTGGACAAGGGAGGGGAAGGGGTTCTTATTACTGACGCGGGTAGCCCCTACTGCTGTGTGGTTCCCCTATTTTTTTTTTTTTCTTTTTGAGACGGAGTCTCGCTCTGTCACCCAGGCTGGAGTGCAGTGGCACAATCTCGGCTCACTGCAAGCTCCACC",
            expected_alt: b"ACG",
        }]);
    }

    // ------------------------------------------------------------------
    // Dot allele / missing: pass-through
    // ------------------------------------------------------------------

    #[test]
    fn test_dot_allele_passthrough() {
        // 20:59 AG -> . : dot allele stays as-is (special)
        // Note: bcftools trims AG/. to A/. but that's because it
        // handles dot specially in multiallelic context. For our
        // single-allele normalize, we skip dot alleles entirely.
        let mut reference = bcftools_norm_reference();
        let result = normalize("20", 59, b"AG", b".", &mut reference).unwrap();
        assert_eq!(result.pos, 59);
        assert_eq!(result.ref_seq, b"AG");
        assert_eq!(result.alt_seq, b".");
        assert!(!result.was_modified);
    }

    // ------------------------------------------------------------------
    // Star allele (spanning deletion): pass-through
    // ------------------------------------------------------------------

    #[test]
    fn test_star_allele_passthrough() {
        // 4:25 T -> * : star allele stays as-is
        let mut reference = bcftools_norm_reference();
        let result = normalize("4", 25, b"T", b"*", &mut reference).unwrap();
        assert_eq!(result.pos, 25);
        assert_eq!(result.ref_seq, b"T");
        assert_eq!(result.alt_seq, b"*");
        assert!(!result.was_modified);
    }

    // ------------------------------------------------------------------
    // SNV multi-allelic (each allele normalized independently)
    // ------------------------------------------------------------------

    #[test]
    fn test_snv_multiallelic_passthrough() {
        // 20:275 A -> C and A -> G: both are SNVs, pass-through
        let mut reference = bcftools_norm_reference();
        let r1 = normalize("20", 275, b"A", b"C", &mut reference).unwrap();
        assert_eq!(r1.pos, 275);
        assert!(!r1.was_modified);

        let r2 = normalize("20", 275, b"A", b"G", &mut reference).unwrap();
        assert_eq!(r2.pos, 275);
        assert!(!r2.was_modified);
    }

    // ------------------------------------------------------------------
    // Indel at position 3 on chr20 (no left-alignment needed)
    // ------------------------------------------------------------------

    #[test]
    fn test_indel_chr20_pos3() {
        // 20:3 G -> CT: already normalized (complex, different lengths)
        run_test_cases(&[NormTestCase {
            name: "indel_chr20_pos3",
            source: "bcftools:norm.vcf:20:3 G/CT (MIT, github.com/samtools/bcftools)",
            chrom: "20",
            pos: 3,
            ref_seq: b"G",
            alt_seq: b"CT",
            expected_pos: 3,
            expected_ref: b"G",
            expected_alt: b"CT",
        }]);
    }

    // ------------------------------------------------------------------
    // Complex multiallelic on chr3 (stays the same)
    // ------------------------------------------------------------------

    #[test]
    fn test_complex_multiallelic_chr3() {
        // 3:10 GTGGAC -> various: each allele individually
        // GTGGAC -> GTGG (deletion): right-trim AC from both, then
        // left-trim GTG from both => pos 13, GAC -> G? No...
        // Let's check: REF=GTGGAC ALT=GTGG
        // Right-trim: no common suffix (C vs G)
        // Left-trim: GTG matches => pos 13, REF=GAC, ALT=G
        // Left-align: indel (ALT len 1), last bases: C vs G, don't match => done
        // So: pos 13, GAC -> G
        //
        // But bcftools output shows 3:10 GTGGAC -> GTGG unchanged.
        // That's because bcftools norm treats the entire multiallelic record
        // as a unit. Our per-allele normalization would trim differently.
        // Let's test per-allele behavior.
        let mut reference = bcftools_norm_reference();

        // GTGGAC -> GTGG: trim left GTG => pos 13, GAC -> G
        let r = normalize("3", 10, b"GTGGAC", b"GTGG", &mut reference).unwrap();
        assert_eq!(r.pos, 13, "GTGGAC/GTGG pos");
        assert_eq!(r.ref_seq, b"GAC", "GTGGAC/GTGG ref");
        assert_eq!(r.alt_seq, b"G", "GTGGAC/GTGG alt");
    }

    // ------------------------------------------------------------------
    // Poly-T region on chr4: deletion with left-alignment
    // ------------------------------------------------------------------

    #[test]
    fn test_poly_t_deletion_chr4() {
        // 4:21 ATTTTTTTTTTTTTTTC -> ATTTTTTTTTTTTTTC (1 T deletion)
        // Right-trim C => ATTTTTTTTTTTTTTT / ATTTTTTTTTTTTTT
        // Left-trim A => pos 22, TTTTTTTTTTTTTTT / TTTTTTTTTTTTTT
        // Left-trim more T's => pos 36, TT / T
        // Left-align: indel (ALT len 1), last bases T==T, prepend upstream
        //   pos 35: ref[34]=T => TT/T... still match => continue
        //   ... all the way back through poly-T run
        // Eventually: pos 22, TT -> T (left edge of poly-T)
        //
        // Actually let's be more careful. After right-trim C:
        // REF = ATTTTTTTTTTTTTTT (16 chars), ALT = ATTTTTTTTTTTTTT (15 chars)
        // Left-trim: A matches => pos 22, REF=TTTTTTTTTTTTTTT (15 T's), ALT=TTTTTTTTTTTTTT (14 T's)
        // Left-trim more T's (while both > 1): trim 13 more T's => pos 35, TT / T
        // Left-align: indel. Last bases T == T. Prepend upstream:
        //   pos 34 (0-based 33): ref="T" => pos 34, REF=TT, ALT=T. Still T==T.
        //   ... continues to pos 22. ref[21] (0-based)="T"
        //   pos 22: REF=TT, ALT=T. Last bases T==T.
        //   pos 21: ref[20] (0-based)="C" => pos 21, REF=CT, ALT=C. C!=T, stop.
        //   Wait, that gives pos 21, CT/C. But there's left-trim: C==C, trim => pos 22, T/"".
        //   Empty ALT isn't valid, so we need the anchor base. Let me re-think.
        //
        // Actually the loop converges across iterations. After first pass:
        // left-align gives us some position with CT/C or similar.
        // Second pass: right-trim nothing, left-trim C => pos+1, T/empty.
        // Empty isn't allowed... the algorithm keeps 1 base minimum.
        //
        // Let me just test what bcftools expects:
        // bcftools output for 4:21 shows it gets decomposed into:
        //   4:25 T/TT,* and 4:36 TC/C,TT,TTC
        // But that's multiallelic split, not what we do per-allele.
        //
        // For the single deletion allele ATTTTTTTTTTTTTTTC -> ATTTTTTTTTTTTTTC:
        // This is deletion of 1 T in a poly-T run. Should left-align to the
        // leftmost position.
        let mut reference = bcftools_norm_reference();

        // Single allele: delete 1 T
        let r = normalize(
            "4",
            21,
            b"ATTTTTTTTTTTTTTTC",
            b"ATTTTTTTTTTTTTTC",
            &mut reference,
        )
        .unwrap();
        // Should left-align through the poly-T run
        // chr4 sequence: TCCCCTCTTGACCTCTCTCTATTTTTTTTTTTTTTTCTGAG...
        //                 1234567890123456789012345678901234567890
        // pos 21 (1-based) = 'A' (0-based index 20)
        // The T-run starts at pos 22 (0-based 21)
        // Left-align should push to pos 21: AT -> A (anchor A before T-run)
        assert_eq!(r.pos, 21, "poly-T del pos");
        assert_eq!(r.ref_seq, b"AT", "poly-T del ref");
        assert_eq!(r.alt_seq, b"A", "poly-T del alt");
    }

    // ------------------------------------------------------------------
    // Insertion of T in poly-T run on chr4
    // ------------------------------------------------------------------

    #[test]
    fn test_poly_t_insertion_chr4() {
        // 4:21 ATTTTTTTTTTTTTTTC -> ATTTTTTTTTTTTTTTTC (1 T insertion)
        // Same poly-T logic, should left-align to pos 21: A -> AT
        let mut reference = bcftools_norm_reference();
        let r = normalize(
            "4",
            21,
            b"ATTTTTTTTTTTTTTTC",
            b"ATTTTTTTTTTTTTTTTC",
            &mut reference,
        )
        .unwrap();
        assert_eq!(r.pos, 21, "poly-T ins pos");
        assert_eq!(r.ref_seq, b"A", "poly-T ins ref");
        assert_eq!(r.alt_seq, b"AT", "poly-T ins alt");
    }

    // ------------------------------------------------------------------
    // Left-align at chromosome boundary (pos 1)
    // ------------------------------------------------------------------

    #[test]
    fn test_left_align_stops_at_pos_1() {
        // If we try to left-align past position 1, it should stop.
        // Create a reference where the variant is at the very start.
        let mut reference = InlineReference::new(HashMap::new()).with_contig("test", b"AAACGT");
        // pos 3: AA -> A (delete 1 A in a poly-A that starts at pos 1)
        let r = normalize("test", 3, b"AA", b"A", &mut reference).unwrap();
        // Should left-align to pos 1: AA -> A
        assert_eq!(r.pos, 1);
        assert_eq!(r.ref_seq, b"AA");
        assert_eq!(r.alt_seq, b"A");
    }

    // ------------------------------------------------------------------
    // Already normalized: was_modified flag
    // ------------------------------------------------------------------

    #[test]
    fn test_was_modified_flag() {
        let mut reference = bcftools_norm_reference();

        // Already normalized SNV
        let r = normalize("20", 81, b"A", b"C", &mut reference).unwrap();
        assert!(!r.was_modified, "SNV should not be modified");

        // Needs normalization
        let r = normalize("20", 80, b"CACAG", b"CACAT", &mut reference).unwrap();
        assert!(r.was_modified, "overspecified SNV should be modified");
    }

    // ------------------------------------------------------------------
    // Symbolic alleles: pass-through
    // ------------------------------------------------------------------

    #[test]
    fn test_symbolic_allele_del() {
        let mut reference = bcftools_norm_reference();
        let r = normalize("20", 24, b"A", b"<DEL>", &mut reference).unwrap();
        assert_eq!(r.pos, 24);
        assert_eq!(r.alt_seq, b"<DEL>");
        assert!(!r.was_modified);
    }

    #[test]
    fn test_symbolic_allele_ins() {
        let mut reference = bcftools_norm_reference();
        let r = normalize("20", 24, b"A", b"<INS>", &mut reference).unwrap();
        assert_eq!(r.pos, 24);
        assert!(!r.was_modified);
    }

    // ------------------------------------------------------------------
    // Right-trim multiallelic (single allele at a time)
    // ------------------------------------------------------------------

    #[test]
    fn test_right_trim_multiallelic_single() {
        // 1:105 TAAACCCTAAA -> TAA: right-trim AA => TAAACCCTA -> T at pos 105
        // Wait, let's be precise:
        // REF = TAAACCCTAAA, ALT = TAA
        // Right-trim: last bases A == A => trim. REF=TAAACCCTA, ALT=TA. A==A => trim.
        // REF=TAAACCCT, ALT=T. T==T => trim. But ALT would become empty (len 0).
        // No: ALT len is 1, so we stop (keep at least 1 base each).
        // Actually "T" has len 1 already, so the while loop condition (a.len() > 1) fails.
        // So after right-trim: REF=TAAACCCTA, ALT=T (trimmed AA from end).
        // Hmm, let's recount. REF=TAAACCCTAAA (11 chars), ALT=TAA (3 chars).
        // Right-trim: A==A (both end in A): TAAACCCTAA / TA (10/2)
        // A==A again: TAAACCCTA / T (9/1). ALT is now len 1, stop.
        // Left-trim: T==T, but ALT is len 1, can't trim (keep ≥1 each).
        // Left-align: indel (ALT len 1). Last bases: A vs T, don't match. Done.
        // Result: pos 105, TAAACCCTA -> T
        run_test_cases(&[NormTestCase {
            name: "right_trim_multiallelic_alt1",
            source: "bcftools:norm.vcf:1:105 TAAACCCTAAA/TAA (MIT, github.com/samtools/bcftools)",
            chrom: "1",
            pos: 105,
            ref_seq: b"TAAACCCTAAA",
            alt_seq: b"TAA",
            expected_pos: 105,
            expected_ref: b"TAAACCCTA",
            expected_alt: b"T",
        }]);
    }

    #[test]
    fn test_right_trim_multiallelic_alt2() {
        // 1:105 TAAACCCTAAA -> TAACCCTAAA: right-trim A => TAAACCCTA / TAACCCTA (10/9)
        // Hmm let me recount. REF=TAAACCCTAAA (11), ALT=TAACCCTAAA (10)
        // Right-trim: A==A => 10/9. A==A => 9/8. A==A => 8/7. Now:
        //   REF=TAAACCCT (8), ALT=TAACCCT (7). T==T => 7/6.
        //   REF=TAAACCC (7), ALT=TAACCC (6). C==C => 6/5.
        //   REF=TAAACC (6), ALT=TAACC (5). C==C => 5/4.
        //   REF=TAAAC (5), ALT=TAAC (4). C==C? No: REF ends in C, ALT ends in C. Yes.
        //   4/3. REF=TAAA (4), ALT=TAA (3). A==A => 3/2. REF=TAA (3), ALT=TA (2).
        //   A==A => 2/1. REF=TA, ALT=T. ALT len 1, stop.
        // Left-trim: T==T but ALT is len 1, stop.
        // Left-align: indel (ALT len 1). Last: A vs T, don't match. Done.
        // Wait that doesn't match bcftools output. Let me re-read norm.out:
        // bcftools shows: 1:105 TAAACCCTA -> T,TAACCCTA
        // So for alt2 (TAACCCTAAA): right-trim AA only (match with REF ending AAA).
        // Wait, the trimming is per-allele in our implementation.
        // REF=TAAACCCTAAA (11 bases), ALT=TAACCCTAAA (10 bases)
        // Right-trim while last bases match:
        //   both end in A: trim => TAAACCCTAA / TAACCCTAA (10/9)
        //   both end in A: trim => TAAACCCTA / TAACCCTA (9/8)
        //   REF ends in A, ALT ends in A: trim => TAAACCCT / TAACCCT (8/7)
        //   T==T: trim => TAAACCC / TAACCC (7/6)
        //   C==C: trim => TAAACC / TAACC (6/5)
        //   C==C: trim => TAAAC / TAAC (5/4)
        //   C==C: NO! REF[4]=C, ALT[3]=C. Wait: TAAAC ends in C, TAAC ends in C. Yes.
        //   trim => TAAA / TAA (4/3)
        //   A==A: trim => TAA / TA (3/2)
        //   A==A: trim => TA / T (2/1). ALT len 1, stop.
        //
        // So right-trim gives: TA / T at pos 105.
        // Left-trim: T==T, but ALT len 1 => stop.
        // Left-align: indel. Last bases A vs T, no match. Done.
        // Result: 105, TA -> T. This is a single-base deletion.
        //
        // But bcftools multiallelic output shows TAACCCTA at pos 105 for alt2.
        // That's because bcftools trims all alleles together to maintain
        // consistency. Our per-allele approach gives different (but equally
        // correct) normalization.
        //
        // Actually wait — per-allele normalization of this is simply
        // "delete 1 A from a run of A's in the REF". The result should
        // left-align to the leftmost A.
        //
        // Let's just check: TA -> T at pos 105. Left-align: ALT is "T",
        // REF is "TA". Last bases: A vs T, don't match. Done.
        // Result: 105, TA -> T.
        //
        // Hmm but actually this should be further left-aligned. The A
        // being deleted is one of many A's. After right-trimming to TA/T,
        // we should left-align: last base of REF=A, last base of ALT=T.
        // They don't match, so no left-alignment. That's correct because
        // the deletion anchor is T (keeps T, deletes A after it).

        run_test_cases(&[NormTestCase {
            name: "right_trim_multiallelic_alt2",
            source:
                "bcftools:norm.vcf:1:105 TAAACCCTAAA/TAACCCTAAA (MIT, github.com/samtools/bcftools)",
            chrom: "1",
            pos: 105,
            ref_seq: b"TAAACCCTAAA",
            alt_seq: b"TAACCCTAAA",
            expected_pos: 105,
            expected_ref: b"TA",
            expected_alt: b"T",
        }]);
    }
}
