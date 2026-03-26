//! Integration tests for variant normalization using the vt normalize
//! test suite (194 indel variants on chr20, b37 assembly).
//!
//! These tests require a real b37 chromosome 20 reference file.  They
//! are skipped at runtime when no reference is available, so they never
//! cause failures in CI or on developer machines that lack the file.
//!
//! To run these tests, place one of the following in your UVID data
//! directory (or set `UVID_DATA_DIR`):
//!
//! - `GRCh37.2bit`  (preferred, ~800 MB)
//! - `GRCh37.fa` + `GRCh37.fa.fai`
//!
//! Source: vt normalize test suite (MIT license, https://github.com/atks/vt)
//! Reference: Adrian Tan et al., Bioinformatics 31(13):2202-2204, 2015
//! https://doi.org/10.1093/bioinformatics/btv112

#[cfg(test)]
mod tests {
    use crate::normalize::algorithm::normalize;
    use crate::normalize::data::open_reference_for_assembly;
    use crate::normalize::reference::ReferenceGenome;

    /// Attempt to open the b37 reference for chr20.
    ///
    /// The vt test data uses bare chromosome names (e.g. "20") which
    /// correspond to the b37/GRCh37 assembly.  UCSC `.2bit` files use
    /// chr-prefixed names ("chr20"), so we use that when querying the
    /// reference.  We try both "GRCh37" and "b37" assembly names.
    ///
    /// Returns `None` if no reference file is found, causing the test
    /// to be skipped rather than failed.
    fn try_open_b37_reference() -> Option<Box<dyn ReferenceGenome>> {
        // Try GRCh37 first (our canonical name), then b37
        open_reference_for_assembly("GRCh37")
            .or_else(|_| open_reference_for_assembly("b37"))
            .ok()
    }

    /// The UCSC .2bit reference uses chr-prefixed contig names regardless
    /// of assembly, so we always query with "chr20" even though the vt
    /// test data uses the bare b37 name "20".
    const VT_CHROM: &str = "chr20";

    /// All 194 test cases from vt normalize test/normalize/01_IN.vcf
    /// and test/normalize/01_OUT.vcf.
    ///
    /// Source: vt (MIT license, https://github.com/atks/vt)
    /// All variants are on chromosome "20" (b37 naming).
    const VT_TEST_CASES: &[(u32, &str, &str, u32, &str, &str)] = &[
        // (in_pos, in_ref, in_alt, exp_pos, exp_ref, exp_alt)
        (421808, "A", "ACCA", 421805, "T", "TCCA"),
        (1292033, "C", "CTTGT", 1292033, "C", "CTTGT"),
        (1340527, "T", "TGTC", 1340527, "T", "TGTC"),
        (1600125, "GAA", "G", 1600125, "GAA", "G"),
        (1728298, "G", "GT", 1728298, "G", "GT"),
        (2171402, "T", "TA", 2171402, "T", "TA"),
        (2171404, "A", "AA", 2171402, "T", "TA"),
        (2982245, "CT", "C", 2982245, "CT", "C"),
        (3025866, "TCAAA", "T", 3025866, "TCAAA", "T"),
        (3373441, "TCTTT", "T", 3373437, "GCTTT", "G"),
        (3635159, "T", "TT", 3635158, "A", "AT"),
        (
            4422119,
            "GCTCCCAGGCTACAGAAAGATGATGGAG",
            "G",
            4422115,
            "GGGAGCTCCCAGGCTACAGAAAGATGAT",
            "G",
        ),
        (5151108, "GTTCT", "G", 5151108, "GTTCT", "G"),
        (5280839, "T", "TATA", 5280839, "T", "TATA"),
        (5291223, "TCAG", "T", 5291223, "TCAG", "T"),
        (5509358, "T", "TG", 5509358, "T", "TG"),
        (5900669, "G", "GC", 5900669, "G", "GC"),
        (5900670, "C", "CC", 5900669, "G", "GC"),
        (6351757, "C", "CTT", 6351757, "C", "CTT"),
        (6362163, "GC", "G", 6362163, "GC", "G"),
        (6481086, "T", "TTGTC", 6481086, "T", "TTGTC"),
        (8080280, "GTTTG", "G", 8080272, "CTTTG", "C"),
        (8781394, "AA", "A", 8781391, "CA", "C"),
        (8833756, "TT", "T", 8833755, "CT", "C"),
        (9035330, "T", "TT", 9035326, "A", "AT"),
        (9311904, "TGTATCTGTCCA", "T", 9311904, "TGTATCTGTCCA", "T"),
        (9389232, "GGGTTTGAT", "G", 9389232, "GGGTTTGAT", "G"),
        (11213808, "A", "AAGCAC", 11213808, "A", "AAGCAC"),
        (11777597, "C", "CT", 11777597, "C", "CT"),
        (12611892, "TCAAGT", "T", 12611888, "TAAGTC", "T"),
        (12656722, "C", "CA", 12656722, "C", "CA"),
        (12664095, "TCA", "T", 12664095, "TCA", "T"),
        (12840732, "A", "ATT", 12840732, "A", "ATT"),
        (12840733, "C", "CT", 12840733, "C", "CT"),
        (13073492, "G", "GT", 13073492, "G", "GT"),
        (13149654, "AATGT", "A", 13149654, "AATGT", "A"),
        (13527241, "T", "TTGT", 13527237, "C", "CTTG"),
        (14563715, "C", "CAT", 14563715, "C", "CAT"),
        (14908890, "T", "TATC", 14908890, "T", "TATC"),
        (14923414, "CTAAAAGCC", "C", 14923410, "GAGCCTAAA", "G"),
        (15082253, "GGTGG", "G", 15082252, "TGGTG", "T"),
        (15317286, "C", "CT", 15317286, "C", "CT"),
        (15524178, "GTTCT", "G", 15524178, "GTTCT", "G"),
        (15701887, "ATCTGAAG", "A", 15701887, "ATCTGAAG", "A"),
        (15701890, "TGAAGTCT", "T", 15701887, "ATCTGAAG", "A"),
        (
            16332171,
            "CAGCTGAGCGCAG",
            "C",
            16332171,
            "CAGCTGAGCGCAG",
            "C",
        ),
        (16336431, "TGAG", "T", 16336431, "TGAG", "T"),
        (16464660, "C", "CCAGT", 16464660, "C", "CCAGT"),
        (16478600, "AAC", "A", 16478600, "AAC", "A"),
        (17187630, "C", "CT", 17187630, "C", "CT"),
        (
            17660514,
            "CTCTCCTAAACCC",
            "C",
            17660506,
            "TCTAAACCCTCTC",
            "T",
        ),
        (18123968, "GAA", "G", 18123968, "GAA", "G"),
        (18487146, "GAATA", "G", 18487146, "GAATA", "G"),
        (18487147, "AATAA", "A", 18487146, "GAATA", "G"),
        (18624883, "TCT", "T", 18624882, "ATC", "A"),
        (18703070, "CT", "C", 18703070, "CT", "C"),
        (19076746, "CTCCCC", "C", 19076743, "ACCCTC", "A"),
        (19250399, "C", "CAT", 19250399, "C", "CAT"),
        (19319754, "T", "TAGCTCCCAACA", 19319754, "T", "TAGCTCCCAACA"),
        (19560555, "GA", "G", 19560555, "GA", "G"),
        (19860444, "GT", "G", 19860444, "GT", "G"),
        (19948252, "C", "CA", 19948252, "C", "CA"),
        (20015197, "A", "AAACA", 20015192, "T", "TAAAC"),
        (20015650, "CTC", "C", 20015647, "ACT", "A"),
        (20033568, "T", "TAATT", 20033563, "A", "ATAAT"),
        (20158176, "TT", "T", 20158175, "CT", "C"),
        (20978417, "T", "TCT", 20978414, "A", "ATC"),
        (21383368, "TACT", "T", 21383366, "TCTA", "T"),
        (21621614, "ATATC", "A", 21621614, "ATATC", "A"),
        (21855321, "C", "CT", 21855321, "C", "CT"),
        (21855322, "T", "TT", 21855321, "C", "CT"),
        (21864933, "GTAATG", "G", 21864931, "CTGTAA", "C"),
        (22088984, "G", "GG", 22088982, "T", "TG"),
        (22414326, "GT", "G", 22414326, "GT", "G"),
        (22449570, "C", "CCTC", 22449567, "T", "TCTC"),
        (22703260, "A", "AGAAG", 22703260, "A", "AGAAG"),
        (22703265, "G", "GAAGG", 22703260, "A", "AGAAG"),
        (22790009, "CTT", "C", 22790009, "CTT", "C"),
        (22984189, "AA", "A", 22984188, "CA", "C"),
        (23635060, "G", "GCGGC", 23635060, "G", "GCGGC"),
        (24183813, "GCTC", "G", 24183813, "GCTC", "G"),
        (24375706, "CAGGATGC", "C", 24375698, "CCAGGATG", "C"),
        (24933608, "C", "CTT", 24933608, "C", "CTT"),
        (
            25179366,
            "GCCACGTCTCTTCTC",
            "G",
            25179366,
            "GCCACGTCTCTTCTC",
            "G",
        ),
        (25258739, "T", "TCTC", 25258739, "T", "TCTC"),
        (25457049, "GCTCCCA", "G", 25457049, "GCTCCCA", "G"),
        (29900012, "GAG", "G", 29900010, "TAG", "T"),
        (30110429, "GTG", "G", 30110427, "CTG", "C"),
        (30158578, "CCA", "C", 30158578, "CCA", "C"),
        (30747542, "CATT", "C", 30747542, "CATT", "C"),
        (30747545, "TATT", "T", 30747542, "CATT", "C"),
        (31031590, "TA", "T", 31031590, "TA", "T"),
        (31403394, "T", "TG", 31403394, "T", "TG"),
        (32537183, "C", "CCT", 32537183, "C", "CCT"),
        (32646307, "A", "ATACT", 32646307, "A", "ATACT"),
        (32871422, "GG", "G", 32871421, "AG", "A"),
        (33626399, "C", "CAGAA", 33626399, "C", "CAGAA"),
        (33626401, "G", "GAAG", 33626399, "C", "CAGA"),
        (35525288, "G", "GGACAG", 35525287, "T", "TGGACA"),
        (35958078, "C", "CTT", 35958078, "C", "CTT"),
        (36010629, "G", "GCAGGGTG", 36010625, "A", "AGGTGCAG"),
        (36521334, "GGTCCAC", "G", 36521334, "GGTCCAC", "G"),
        (36686810, "CTT", "C", 36686810, "CTT", "C"),
        (36686811, "TTT", "T", 36686810, "CTT", "C"),
        (36791656, "A", "ACT", 36791656, "A", "ACT"),
        (37037631, "AT", "A", 37037631, "AT", "A"),
        (37048678, "A", "ATA", 37048675, "A", "AAT"),
        (37388010, "ACTCTAGCGAGA", "A", 37388009, "AACTCTAGCGAG", "A"),
        (37394795, "AT", "A", 37394795, "AT", "A"),
        (37394796, "TT", "T", 37394795, "AT", "A"),
        (37706107, "A", "AC", 37706107, "A", "AC"),
        (38022703, "GAAG", "G", 38022700, "CAAG", "C"),
        (38258947, "G", "GATGGA", 38258947, "G", "GATGGA"),
        (38416239, "TTGAT", "T", 38416231, "CTGAT", "C"),
        (38441280, "AA", "A", 38441276, "TA", "T"),
        (39074861, "TCT", "T", 39074860, "TTC", "T"),
        (39475162, "TG", "T", 39475162, "TG", "T"),
        (39611561, "AG", "A", 39611561, "AG", "A"),
        (39918080, "AGA", "A", 39918077, "CAG", "C"),
        (40552038, "TATG", "T", 40552038, "TATG", "T"),
        (41086660, "AA", "A", 41086658, "CA", "C"),
        (41284374, "G", "GTT", 41284374, "G", "GTT"),
        (41284377, "T", "TTT", 41284374, "G", "GTT"),
        (41771255, "CAA", "C", 41771255, "CAA", "C"),
        (41900315, "GT", "G", 41900315, "GT", "G"),
        (42493962, "T", "TA", 42493962, "T", "TA"),
        (42551667, "C", "CAG", 42551667, "C", "CAG"),
        (42551669, "G", "GAG", 42551667, "C", "CAG"),
        (42595607, "ACCAACA", "A", 42595599, "CCACCAA", "C"),
        (43521315, "ACTCTA", "A", 43521313, "TTACTC", "T"),
        (43893147, "A", "AGT", 43893147, "A", "AGT"),
        (44442298, "C", "CTA", 44442298, "C", "CTA"),
        (45088514, "AGATGATATT", "A", 45088514, "AGATGATATT", "A"),
        (45849248, "T", "TAACA", 45849248, "T", "TAACA"),
        (45850737, "GAG", "G", 45850735, "TAG", "T"),
        (45924828, "TT", "T", 45924827, "CT", "C"),
        (46657104, "C", "CTCC", 46657103, "A", "ACTC"),
        (46687122, "T", "TTA", 46687122, "T", "TTA"),
        (46687123, "T", "TAT", 46687122, "T", "TTA"),
        (46981903, "GATGAGCTCTAA", "G", 46981903, "GATGAGCTCTAA", "G"),
        (46981904, "ATGAGCTCTAAA", "A", 46981903, "GATGAGCTCTAA", "G"),
        (47673119, "CCAA", "C", 47673119, "CCAA", "C"),
        (48602568, "AG", "A", 48602568, "AG", "A"),
        (48696636, "G", "GG", 48696635, "A", "AG"),
        (49350131, "A", "ATA", 49350127, "C", "CTA"),
        (
            50018388,
            "GGTGTCAACAAATAG",
            "G",
            50018386,
            "AAGGTGTCAACAAAT",
            "A",
        ),
        (50114687, "TT", "T", 50114685, "CT", "C"),
        (50711979, "TG", "T", 50711979, "TG", "T"),
        (50712021, "C", "CTC", 50712019, "G", "GTC"),
        (50897361, "G", "GGAATGTCAGCC", 50897361, "G", "GGAATGTCAGCC"),
        (51027772, "CCT", "C", 51027772, "CCT", "C"),
        (51055421, "AC", "A", 51055421, "AC", "A"),
        (51127873, "ATTCCAAATCC", "A", 51127873, "ATTCCAAATCC", "A"),
        (51235290, "C", "CTT", 51235290, "C", "CTT"),
        (51512776, "TCT", "T", 51512772, "GCT", "G"),
        (51630299, "AC", "A", 51630299, "AC", "A"),
        (51731206, "C", "CC", 51731202, "T", "TC"),
        (52772453, "A", "AA", 52772452, "T", "TA"),
        (53713059, "AG", "A", 53713059, "AG", "A"),
        (53802839, "TATAG", "T", 53802839, "TATAG", "T"),
        (53804546, "AAC", "A", 53804546, "AAC", "A"),
        (54125169, "T", "TAA", 54125169, "T", "TAA"),
        (54158733, "T", "TCTCT", 54158730, "G", "GTCTC"),
        (54375029, "A", "ACACT", 54375029, "A", "ACACT"),
        (54664038, "A", "AAA", 54664037, "T", "TAA"),
        (54665425, "T", "TTATTTCCA", 54665425, "T", "TTATTTCCA"),
        (55292357, "GTT", "G", 55292357, "GTT", "G"),
        (55292358, "TTT", "T", 55292357, "GTT", "G"),
        (55444416, "GGAG", "G", 55444413, "AGAG", "A"),
        (55557713, "AGCTCTGA", "A", 55557711, "AGAGCTCT", "A"),
        (55622401, "CAG", "C", 55622401, "CAG", "C"),
        (56311463, "T", "TT", 56311461, "A", "AT"),
        (56383576, "TTTG", "T", 56383576, "TTTG", "T"),
        (56564760, "C", "CT", 56564760, "C", "CT"),
        (56788029, "T", "TTA", 56788029, "T", "TTA"),
        (56788032, "T", "TAT", 56788029, "T", "TTA"),
        (57863967, "A", "AAACCA", 57863961, "G", "GAAACC"),
        (
            57882335,
            "T",
            "TCATGATGCCAAG",
            57882335,
            "T",
            "TCATGATGCCAAG",
        ),
        (57967774, "T", "TT", 57967772, "A", "AT"),
        (58475138, "AC", "A", 58475138, "AC", "A"),
        (58742718, "TACTGAT", "T", 58742717, "CTACTGA", "C"),
        (59201326, "G", "GT", 59201326, "G", "GT"),
        (
            59641522,
            "C",
            "CAGGGAACCAAGCGAGGGAGATTCAGACCCTGCCT",
            59641522,
            "C",
            "CAGGGAACCAAGCGAGGGAGATTCAGACCCTGCCT",
        ),
        (59705557, "TTTGT", "T", 59705552, "CTTTG", "C"),
        (59917975, "A", "ATTA", 59917971, "C", "CATT"),
        (59961560, "TTT", "T", 59961559, "CTT", "C"),
        (60159925, "TGTC", "T", 60159925, "TGTC", "T"),
        (60274581, "TC", "T", 60274581, "TC", "T"),
        (60378250, "C", "CCA", 60378250, "C", "CCA"),
        (60744904, "ACCGTCCACA", "A", 60744897, "TGTCCACACC", "T"),
        (62523744, "C", "CTT", 62523744, "C", "CTT"),
        (62762721, "C", "CAT", 62762721, "C", "CAT"),
        (62762759, "CAT", "C", 62762759, "CAT", "C"),
        (62812114, "AGC", "A", 62812114, "AGC", "A"),
    ];

    /// Run all 194 vt normalization test cases.
    ///
    /// The test is skipped (not failed) when no b37 reference genome is
    /// available in the data directory.
    #[test]
    fn test_vt_normalize_194_cases() {
        let mut reference = match try_open_b37_reference() {
            Some(r) => r,
            None => {
                eprintln!(
                    "SKIPPED: vt normalization tests require a GRCh37 reference file.\n\
                     Place GRCh37.2bit or GRCh37.fa (+.fai) in the UVID data directory,\n\
                     or set UVID_DATA_DIR to a directory containing the file."
                );
                return;
            }
        };

        let total = VT_TEST_CASES.len();
        let mut passed = 0;
        let mut failed = 0;
        let mut failures: Vec<String> = Vec::new();

        for (i, &(in_pos, in_ref, in_alt, exp_pos, exp_ref, exp_alt)) in
            VT_TEST_CASES.iter().enumerate()
        {
            let result = match normalize(
                VT_CHROM,
                in_pos,
                in_ref.as_bytes(),
                in_alt.as_bytes(),
                reference.as_mut(),
            ) {
                Ok(r) => r,
                Err(e) => {
                    failures.push(format!(
                        "  [{}/{}] {}:{} {}/{} => ERROR: {}",
                        i + 1,
                        total,
                        VT_CHROM,
                        in_pos,
                        in_ref,
                        in_alt,
                        e
                    ));
                    failed += 1;
                    continue;
                }
            };

            let pos_ok = result.pos == exp_pos;
            let ref_ok = result.ref_seq == exp_ref.as_bytes();
            let alt_ok = result.alt_seq == exp_alt.as_bytes();

            if pos_ok && ref_ok && alt_ok {
                passed += 1;
            } else {
                let msg = format!(
                    "  [{}/{}] {}:{} {}/{} => got {}:{}/{}, expected {}:{}/{}",
                    i + 1,
                    total,
                    VT_CHROM,
                    in_pos,
                    in_ref,
                    in_alt,
                    result.pos,
                    String::from_utf8_lossy(&result.ref_seq),
                    String::from_utf8_lossy(&result.alt_seq),
                    exp_pos,
                    exp_ref,
                    exp_alt,
                );
                failures.push(msg);
                failed += 1;
            }
        }

        if !failures.is_empty() {
            panic!(
                "vt normalization: {}/{} passed, {} failed:\n{}",
                passed,
                total,
                failed,
                failures.join("\n")
            );
        }

        eprintln!("vt normalization: all {}/{} cases passed", passed, total);
    }

    /// Verify that the count of test cases that should change matches
    /// the vt expected output (80 out of 194 variants are modified).
    #[test]
    fn test_vt_normalize_modification_count() {
        let mut reference = match try_open_b37_reference() {
            Some(r) => r,
            None => {
                eprintln!("SKIPPED: requires GRCh37 reference");
                return;
            }
        };

        let mut modified_count = 0;
        for &(in_pos, in_ref, in_alt, _, _, _) in VT_TEST_CASES.iter() {
            let result = normalize(
                VT_CHROM,
                in_pos,
                in_ref.as_bytes(),
                in_alt.as_bytes(),
                reference.as_mut(),
            )
            .unwrap_or_else(|e| {
                panic!(
                    "normalize {}:{} {}/{}: {}",
                    VT_CHROM, in_pos, in_ref, in_alt, e
                )
            });
            if result.was_modified {
                modified_count += 1;
            }
        }

        // vt reports 80 variants modified out of 194
        assert_eq!(
            modified_count, 80,
            "expected 80 modified variants (matching vt output), got {}",
            modified_count
        );
    }
}
