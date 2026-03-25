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
//! Source: GATK LeftAlignAndTrimVariants test suite
//!   (Apache License 2.0, github.com/broadinstitute/gatk)
//! Reference: Broad Institute, Genome Analysis Toolkit (GATK)

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
    /// Source: GATK (Apache License 2.0, github.com/broadinstitute/gatk)
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

        // Case 4a: Multiallelic decomposed — insertion allele (already normalized)
        ("chr20", 19883345, "T", "TT",
            19883345, "T", "TT",
            "multiallelic insertion allele, already normalized"),

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

        // Case 9: 296bp insertion, already left-aligned
        ("chr20", 63669973, "G",
            "GGACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG",
            63669973, "G",
            "GGACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG\
             GACAGACGTTTCGCCAAGATGGGTGGAATGGCCAGTTAACCACTGGGAGAGCATCCG",
            "296bp insertion already left-aligned"),

        // Case 10: 296bp deletion, already left-aligned
        ("chr20", 64012187,
            "TACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT",
            "T",
            64012187,
            "TACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT\
             ACACCTACGAGAGGAGGACAGGAAGTCTCACGGCTTCTTTCTTTGCTAGTAGAAATGAT",
            "T",
            "296bp deletion already left-aligned"),

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

        // Case 16: 2bp deletion that left-aligns past a nearby SNP (tricky!)
        //
        // GATK notes this as "particularly tricky": the deletion AAA→A at
        // chr21:13255301 left-aligns to AAA→A at chr21:13255297, passing a
        // SNP A→G at chr21:13255296.
        ("chr21", 13255301, "AAA", "A",
            13255297, "AAA", "A",
            "deletion left-aligns past SNP (tricky)"),

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
    /// - Already normalized: 5 unchanged (cases 4a, 7a, 9, 10, 14)
    /// - Modified: 10 (cases 1, 2, 5, 6, 8, 11, 12, 13, 16, 17)
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

        // 10 of 19 per-allele cases should be modified
        assert_eq!(
            modified_count, 10,
            "expected 10 modified alleles (GATK left-align+trim), got {}",
            modified_count
        );
    }
}
