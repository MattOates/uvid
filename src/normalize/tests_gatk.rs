//! Integration tests for variant normalization using the GATK
//! LeftAlignAndTrimVariants test suite (17 variants on chr20/chr21,
//! GRCh38 assembly).
//!
//! These tests require a real GRCh38 reference file.  They are skipped
//! at runtime when no reference is available, so they never cause
//! failures in CI or on developer machines that lack the file.
//!
//! To run these tests, place one of the following in your UVID data
//! directory (or set `UVID_DATA_DIR`):
//!
//! - `GRCh38.2bit`  (preferred, ~800 MB)
//! - `GRCh38.fa` + `GRCh38.fa.fai`
//!
//! ## Divergence from GATK expected output
//!
//! Our normalizer implements the Tan et al. 2015 algorithm (the same
//! approach used by bcftools and vt), which fully left-aligns every
//! indel to the leftmost possible position.  GATK's
//! `LeftAlignAndTrimVariants` has two known limitations that cause it
//! to stop short:
//!
//! 1. **`--max-indel-length` (default 200)** — GATK skips
//!    normalization entirely for indels whose allele length exceeds
//!    this threshold.  Our normalizer has no such limit.
//!    Affects GATK test cases 9 and 10 (296 bp insertion/deletion).
//!
//! 2. **`distanceToLastVariant` constraint** — GATK refuses to
//!    left-align a variant past the end position of the previously
//!    emitted variant on the same contig.  This prevents full
//!    left-alignment when a decomposed multiallelic or nearby SNP is
//!    present.  The GATK team consider this a bug; see
//!    <https://github.com/broadinstitute/gatk/pull/9305> (Dec 2025,
//!    open) and the kachulis review comments on
//!    <https://github.com/broadinstitute/gatk/pull/6427>.
//!    Affects GATK test cases 4a and 16.
//!
//! In each divergent case the expected values below reflect the
//! *fully* left-aligned result, verified against the GRCh38 reference
//! genome.
//!
//! Source: GATK LeftAlignAndTrimVariants test suite
//!   (Apache License 2.0, https://github.com/broadinstitute/gatk)
//! Reference: Broad Institute, Genome Analysis Toolkit (GATK)
//!   https://gatk.broadinstitute.org/

#[cfg(test)]
mod tests {
    use crate::normalize::algorithm::normalize;
    use crate::normalize::data::open_reference_for_assembly;
    use crate::normalize::reference::ReferenceGenome;

    /// Attempt to open the GRCh38 reference.
    ///
    /// The GATK test data uses chr-prefixed chromosome names (e.g.
    /// "chr20", "chr21") corresponding to the GRCh38 assembly.  We
    /// also try "hg38" as an alternate name.
    ///
    /// Returns `None` if no reference file is found, causing the test
    /// to be skipped rather than failed.
    fn try_open_grch38_reference() -> Option<Box<dyn ReferenceGenome>> {
        open_reference_for_assembly("GRCh38")
            .or_else(|_| open_reference_for_assembly("hg38"))
            .ok()
    }

    /// All single-allele test cases from the GATK LeftAlignAndTrimVariants
    /// test suite (test_left_align_hg38.vcf → expected_left_align_hg38.vcf).
    ///
    /// Multiallelic cases (input cases 4 and 7: T→TT,C and T→TC,C) are
    /// decomposed into individual alleles here, since our implementation
    /// normalizes each allele independently.  Both alleles in those
    /// records are already in normalized form, so expected == input.
    ///
    /// Source: GATK (Apache License 2.0, https://github.com/broadinstitute/gatk)
    ///
    /// Format: (chrom, in_pos, in_ref, in_alt, exp_pos, exp_ref, exp_alt, description)
    #[rustfmt::skip]
    const GATK_TEST_CASES: &[(&str, u32, &str, &str, u32, &str, &str, &str)] = &[
        // Case 1: Large insertion (62bp), left-aligns 132 positions
        ("chr20", 1259782, "G",
            "GATCTTCCCTCTTTTCTAATATAAACACATAAAGCTCTGTTTCCTTCTAGGTAACTGGTTTGAG",
            1259650, "A",
            "ATTTGAGATCTTCCCTCTTTTCTAATATAAACACATAAAGCTCTGTTTCCTTCTAGGTAACTGG",
            "large insertion left-aligns 132bp"),

        // Case 2: 3bp deletion, left-aligns 10 positions
        ("chr20", 19285487, "AAAA", "A",
            19285477, "CAAA", "C",
            "deletion left-aligns 10bp"),

        // Case 3: SNV, pass-through
        ("chr20", 19285500, "A", "C",
            19285500, "A", "C",
            "SNV pass-through"),

        // Case 4a: Multiallelic decomposed — insertion allele
        //
        // GATK divergence: GATK reports this as "already normalized"
        // (T/TT at 19883345) because its `distanceToLastVariant`
        // constraint prevents left-aligning past the previously
        // emitted variant at this site.  Full left-alignment places
        // the T insertion at the start of a poly-T run
        // (19883338–19883359), anchored by the C at 19883337.
        // See: https://github.com/broadinstitute/gatk/pull/9305
        ("chr20", 19883345, "T", "TT",
            19883337, "C", "CT",
            "multiallelic insertion allele, left-aligns in poly-T (GATK stops short)"),

        // Case 4b: Multiallelic decomposed — SNV allele (already normalized)
        ("chr20", 19883345, "T", "C",
            19883345, "T", "C",
            "multiallelic SNV allele, already normalized"),

        // Case 5: Overspecified deletion, trims shared suffix
        ("chr20", 19883389, "GAGT", "GA",
            19883390, "AGT", "A",
            "overspecified deletion trims suffix"),

        // Case 6: Overspecified complex → deletion after suffix trim
        ("chr20", 19883397, "CGGA", "CA",
            19883397, "CGG", "C",
            "overspecified complex trims to deletion"),

        // Case 7a: Multiallelic decomposed — insertion allele (already normalized)
        ("chr20", 19883412, "T", "TC",
            19883412, "T", "TC",
            "multiallelic insertion allele, already normalized"),

        // Case 7b: Multiallelic decomposed — SNV allele (already normalized)
        ("chr20", 19883412, "T", "C",
            19883412, "T", "C",
            "multiallelic SNV allele, already normalized"),

        // Case 8: 5bp deletion, left-aligns 9 positions + trims
        ("chr20", 19885711, "AAGAAAA", "AA",
            19885702, "CGAAAA", "C",
            "deletion left-aligns 9bp and trims"),

        // Case 9: 296bp insertion, left-aligns through tandem repeat
        //
        // GATK divergence: GATK skips this variant entirely because
        // it exceeds `--max-indel-length` (default 200 bp).  Full
        // left-alignment shifts the insertion 329 bp upstream through
        // ~6 copies of a 57 bp tandem repeat, anchored by the A at
        // 63669644.
        ("chr20", 63669973, "G",
            "GGACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG",
            63669644, "A",
            "ACCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCG",
            "296bp insertion left-aligns through tandem repeat (GATK skips >200bp)"),

        // Case 10: 296bp deletion, left-aligns through tandem repeat
        //
        // GATK divergence: same `--max-indel-length` skip as case 9.
        // Full left-alignment shifts the deletion 382 bp upstream
        // through ~6 copies of a ~59 bp tandem repeat, anchored by
        // the T at 64011805.
        ("chr20", 64012187,
            "TACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT",
            "T",
            64011805,
            "TGGCTTCTTTCTTTGCTAGTAGAAATGATACACCTACGAGAGGAGGACAGGAAGTCTCAC\
             GGCTTCTTTCTTTGCTAGTAGAAATGATACACCTACGAGAGGAGGACAGGAAGTCTCAC\
             GGCTTCTTTCTTTGCTAGTAGAAATGATACACCTACGAGAGGAGGACAGGAAGTCTCAC\
             GGCTTCTTTCTTTGCTAGTAGAAATGATACACCTACGAGAGGAGGACAGGAAGTCTCAC\
             GGCTTCTTTCTTTGCTAGTAGAAATGATACACCTACGAGAGGAGGACAGGAAGTCTCAC",
            "T",
            "296bp deletion left-aligns through tandem repeat (GATK skips >200bp)"),

        // Case 11: 6bp insertion in TG repeat, left-aligns 6 positions
        ("chr21", 8405579, "G", "GTGTGTG",
            8405573, "A", "ATGTGTG",
            "insertion in TG repeat left-aligns 6bp"),

        // Case 12: 2bp insertion in poly-T, left-aligns 6 positions
        ("chr21", 10382395, "T", "TTT",
            10382389, "A", "ATT",
            "insertion in poly-T left-aligns 6bp"),

        // Case 13: 3bp deletion in repeat, left-aligns 16 positions
        ("chr21", 10388249, "GAAG", "G",
            10388233, "GGAA", "G",
            "trinucleotide deletion left-aligns 16bp"),

        // Case 14: 2bp insertion, already left-aligned
        ("chr21", 10804284, "T", "TGC",
            10804284, "T", "TGC",
            "GC insertion already left-aligned"),

        // Case 15: SNV, pass-through
        ("chr21", 13255296, "A", "G",
            13255296, "A", "G",
            "SNV pass-through"),

        // Case 16: 2bp deletion that left-aligns past a nearby SNP
        //
        // GATK divergence: GATK stops left-aligning at chr21:13255297
        // because of its `distanceToLastVariant` constraint — the
        // previous variant (SNP A→G at chr21:13255296, case 15) blocks
        // further movement.  Full left-alignment places the deletion
        // at the start of the poly-A run (13255290–13255304), anchored
        // by the C at 13255289.
        // See: https://github.com/broadinstitute/gatk/pull/9305
        //      https://github.com/broadinstitute/gatk/pull/6427
        ("chr21", 13255301, "AAA", "A",
            13255289, "CAA", "C",
            "deletion left-aligns past SNP to poly-A start (GATK stops at SNP)"),

        // Case 17: 44bp deletion in TTCCC repeat, left-aligns 189 positions
        ("chr21", 39584006,
            "CTTCCCTTCCCTTCCCTTCCCTTCCCTTCCCTTCCCTTCCCTTCCC", "C",
            39583817,
            "CTCCCTTCCCTTCCCTTCCCTTCCCTTCCCTTCCCTTCCCTTCCCT", "C",
            "large repeat deletion left-aligns 189bp"),
    ];

    /// Run all GATK test cases against our normalizer.
    ///
    /// Each case is normalized independently.  Cases that require a
    /// reference genome are skipped (with a message) when the GRCh38
    /// reference is not available.
    #[test]
    fn test_gatk_normalize_17_cases() {
        let mut reference = match try_open_grch38_reference() {
            Some(r) => r,
            None => {
                eprintln!("SKIPPED: requires GRCh38 reference");
                return;
            }
        };

        let total = GATK_TEST_CASES.len();
        let mut passed = 0;
        let mut failed = 0;
        let mut failures = Vec::new();

        for (i, &(chrom, in_pos, in_ref, in_alt, exp_pos, exp_ref, exp_alt, desc)) in
            GATK_TEST_CASES.iter().enumerate()
        {
            let result = normalize(
                chrom,
                in_pos,
                in_ref.as_bytes(),
                in_alt.as_bytes(),
                reference.as_mut(),
            )
            .unwrap_or_else(|e| {
                panic!(
                    "case {}: normalize {}:{} {}/{}: {}",
                    i + 1,
                    chrom,
                    in_pos,
                    in_ref,
                    in_alt,
                    e
                )
            });

            let got_ref = std::str::from_utf8(&result.ref_seq).unwrap();
            let got_alt = std::str::from_utf8(&result.alt_seq).unwrap();

            if result.pos != exp_pos || got_ref != exp_ref || got_alt != exp_alt {
                failed += 1;
                failures.push(format!(
                    "  case {} ({}): {}:{} {}/{}  =>  got {}:{} {}/{}  expected {}:{} {}/{}",
                    i + 1,
                    desc,
                    chrom,
                    in_pos,
                    in_ref,
                    in_alt,
                    chrom,
                    result.pos,
                    got_ref,
                    got_alt,
                    chrom,
                    exp_pos,
                    exp_ref,
                    exp_alt,
                ));
            } else {
                passed += 1;
            }
        }

        if failed > 0 {
            panic!(
                "GATK normalization: {}/{} cases FAILED:\n{}",
                failed,
                total,
                failures.join("\n")
            );
        }

        eprintln!("GATK normalization: all {}/{} cases passed", passed, total);
    }

    /// Verify the expected modification count.
    ///
    /// Of the 19 per-allele cases (17 GATK records with 2 multiallelic
    /// records decomposed into 2 alleles each = 19 allele-level cases):
    /// - SNVs: 4 unchanged (cases 3, 4b, 7b, 15)
    /// - Already normalized: 2 unchanged (cases 7a, 14)
    /// - Modified: 13 (cases 1, 2, 4a, 5, 6, 8, 9, 10, 11, 12, 13, 16, 17)
    ///
    /// Note: GATK reports only 10 modified because it skips indels
    /// >200 bp (cases 9, 10) and constrains left-alignment by
    /// `distanceToLastVariant` (case 4a).
    #[test]
    fn test_gatk_normalize_modification_count() {
        let mut reference = match try_open_grch38_reference() {
            Some(r) => r,
            None => {
                eprintln!("SKIPPED: requires GRCh38 reference");
                return;
            }
        };

        let mut modified_count = 0;
        for &(chrom, in_pos, in_ref, in_alt, _, _, _, _) in GATK_TEST_CASES.iter() {
            let result = normalize(
                chrom,
                in_pos,
                in_ref.as_bytes(),
                in_alt.as_bytes(),
                reference.as_mut(),
            )
            .unwrap_or_else(|e| {
                panic!(
                    "normalize {}:{} {}/{}: {}",
                    chrom, in_pos, in_ref, in_alt, e
                )
            });
            if result.was_modified {
                modified_count += 1;
            }
        }

        // 13 of 19 per-allele cases should be modified
        assert_eq!(
            modified_count, 13,
            "expected 13 modified alleles (full left-align+trim), got {}",
            modified_count
        );
    }
}
