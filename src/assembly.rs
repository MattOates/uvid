/// Genome assembly definitions for GRCh37 and GRCh38.
///
/// Provides chromosome metadata and encoding/decoding between
/// chromosome name + position and the linearized genome position
/// used in the UVID 128-bit layout.
///
/// Linearized positions concatenate all chromosomes into a single
/// coordinate space: chr1 starts at offset 0, chr2 at the end of chr1,
/// etc. The linearized position for a variant is `offset[chr] + pos`.
use std::fmt;

/// Supported genome assemblies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Assembly {
    GRCh37 = 0,
    GRCh38 = 1,
}

impl Assembly {
    /// Decode from the 2-bit assembly field in a UVID.
    pub fn from_bits(bits: u8) -> Option<Self> {
        match bits {
            0 => Some(Assembly::GRCh37),
            1 => Some(Assembly::GRCh38),
            _ => None,
        }
    }

    /// Encode to the 2-bit assembly field.
    pub fn to_bits(self) -> u8 {
        self as u8
    }
}

impl fmt::Display for Assembly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Assembly::GRCh37 => write!(f, "GRCh37"),
            Assembly::GRCh38 => write!(f, "GRCh38"),
        }
    }
}

impl std::str::FromStr for Assembly {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "grch37" | "hg19" | "37" => Ok(Assembly::GRCh37),
            "grch38" | "hg38" | "38" => Ok(Assembly::GRCh38),
            _ => Err(format!("Unknown assembly: {}", s)),
        }
    }
}

/// Metadata for a single chromosome within an assembly.
#[derive(Debug, Clone)]
pub struct ChromosomeInfo {
    pub name: &'static str,
    pub length: u32,
    pub genbank: Option<&'static str>,
    pub refseq: &'static str,
}

/// Contig naming convention used in a VCF file's CHROM column.
///
/// VCF files use different naming schemes depending on the reference genome
/// and pipeline that produced them.  Detecting the scheme from the VCF header
/// (or first data record) allows [`ChrIndex::resolve`] to map any CHROM value
/// to a chromosome index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContigScheme {
    /// Bare numeric/letter names: `"1"`, `"2"`, ..., `"X"`, `"Y"`, `"MT"`.
    ///
    /// Used by b37, humanG1Kv37, 1000 Genomes, Ensembl, and GRC assemblies.
    Bare,
    /// UCSC chr-prefixed names: `"chr1"`, `"chr2"`, ..., `"chrX"`, `"chrY"`, `"chrM"`.
    ///
    /// Used by hg19, hg38, Dragen, DeepVariant, and most clinical pipelines.
    Ucsc,
    /// NCBI RefSeq accessions: `"NC_000001.11"`, `"NC_000023.11"`, `"NC_012920.1"`.
    ///
    /// Used by dbSNP VCF, ClinVar VCF, and NCBI RefSeq reference FASTAs.
    /// Version numbers are assembly-specific (e.g. `.11` for GRCh38, `.10`
    /// for GRCh37 on chr1).
    RefSeq,
    /// GenBank accessions: `"CM000663.2"`, `"CM000664.2"`, etc.
    ///
    /// Rare in practice; appears when GenBank assembly FASTAs (`GCA_*`) are
    /// used as the reference.
    GenBank,
}

/// Classify a single CHROM string into its naming scheme.
///
/// Uses prefix-based heuristics (case-insensitive):
/// - `NC_` → [`ContigScheme::RefSeq`]
/// - `CM` followed by a digit → [`ContigScheme::GenBank`]
/// - `chr` → [`ContigScheme::Ucsc`]
/// - anything else → [`ContigScheme::Bare`]
pub fn classify_chrom(name: &str) -> ContigScheme {
    let upper = name.to_ascii_uppercase();
    if upper.starts_with("NC_") {
        ContigScheme::RefSeq
    } else if upper.starts_with("CM") && upper.as_bytes().get(2).is_some_and(|b| b.is_ascii_digit())
    {
        ContigScheme::GenBank
    } else if upper.starts_with("CHR") {
        ContigScheme::Ucsc
    } else {
        ContigScheme::Bare
    }
}

/// Detect the contig naming scheme from a set of VCF contig IDs.
///
/// Classifies each ID and returns the most specific scheme observed.
/// Priority: RefSeq > GenBank > Ucsc > Bare.  If `contig_ids` is empty,
/// returns [`ContigScheme::Bare`] as the default.
pub fn detect_contig_scheme(contig_ids: &[&str]) -> ContigScheme {
    let mut best = ContigScheme::Bare;
    for id in contig_ids {
        let scheme = classify_chrom(id);
        // Return immediately on the most specific schemes
        match scheme {
            ContigScheme::RefSeq => return ContigScheme::RefSeq,
            ContigScheme::GenBank => return ContigScheme::GenBank,
            ContigScheme::Ucsc => best = ContigScheme::Ucsc,
            ContigScheme::Bare => {}
        }
    }
    best
}

/// Chromosome index in the UVID encoding (5-bit field, 0-24).
///
/// chr1=0, chr2=1, ..., chr22=21, chrX=22, chrY=23, chrM=24
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChrIndex(pub u8);

impl ChrIndex {
    /// Maximum valid chromosome index.
    pub const MAX: u8 = 24;

    /// Parse a chromosome name string to its index.
    ///
    /// Accepts: "1"-"22", "X", "Y", "M", "MT",
    /// and "chr1"-"chr22", "chrX", "chrY", "chrM", "chrMT"
    pub fn from_name(name: &str) -> Option<Self> {
        let name = name.strip_prefix("chr").unwrap_or(name);
        match name {
            "1" => Some(ChrIndex(0)),
            "2" => Some(ChrIndex(1)),
            "3" => Some(ChrIndex(2)),
            "4" => Some(ChrIndex(3)),
            "5" => Some(ChrIndex(4)),
            "6" => Some(ChrIndex(5)),
            "7" => Some(ChrIndex(6)),
            "8" => Some(ChrIndex(7)),
            "9" => Some(ChrIndex(8)),
            "10" => Some(ChrIndex(9)),
            "11" => Some(ChrIndex(10)),
            "12" => Some(ChrIndex(11)),
            "13" => Some(ChrIndex(12)),
            "14" => Some(ChrIndex(13)),
            "15" => Some(ChrIndex(14)),
            "16" => Some(ChrIndex(15)),
            "17" => Some(ChrIndex(16)),
            "18" => Some(ChrIndex(17)),
            "19" => Some(ChrIndex(18)),
            "20" => Some(ChrIndex(19)),
            "21" => Some(ChrIndex(20)),
            "22" => Some(ChrIndex(21)),
            "X" | "x" => Some(ChrIndex(22)),
            "Y" | "y" => Some(ChrIndex(23)),
            "M" | "m" | "MT" | "mt" => Some(ChrIndex(24)),
            _ => None,
        }
    }

    /// Convert a chromosome index back to its canonical name (without "chr" prefix).
    pub fn to_name(self) -> Option<&'static str> {
        CHROMOSOME_NAMES.get(self.0 as usize).copied()
    }

    /// Convert a chromosome index to its UCSC-style name (`"chr1"`, `"chrX"`, `"chrM"`).
    ///
    /// This is the naming convention used by UCSC `.2bit` reference genome files
    /// (e.g. `hg38.2bit`).  Use this when querying a reference genome for
    /// normalization.
    pub fn to_ucsc_name(self) -> Option<&'static str> {
        UCSC_CHROMOSOME_NAMES.get(self.0 as usize).copied()
    }

    /// Parse a RefSeq accession to a chromosome index using [`ChromosomeInfo`] data.
    ///
    /// Performs **strict version matching**: the full accession including version
    /// suffix must match exactly (e.g. `"NC_000001.11"` matches chr1 in GRCh38
    /// but not in GRCh37, which uses `"NC_000001.10"`).
    pub fn from_refseq(accession: &str, assembly: Assembly) -> Option<Self> {
        for (i, info) in chromosomes(assembly).iter().enumerate() {
            if info.refseq == accession {
                return Some(ChrIndex(i as u8));
            }
        }
        None
    }

    /// Parse a GenBank accession to a chromosome index using [`ChromosomeInfo`] data.
    ///
    /// Performs **strict version matching**: the full accession including version
    /// suffix must match exactly (e.g. `"CM000663.2"` matches chr1 in GRCh38
    /// but not in GRCh37, which uses `"CM000663.1"`).
    ///
    /// Note: chrM has no GenBank accession (`genbank: None` in [`ChromosomeInfo`]).
    pub fn from_genbank(accession: &str, assembly: Assembly) -> Option<Self> {
        for (i, info) in chromosomes(assembly).iter().enumerate() {
            if info.genbank == Some(accession) {
                return Some(ChrIndex(i as u8));
            }
        }
        None
    }

    /// Resolve a CHROM string to a chromosome index using the detected
    /// naming scheme and assembly context.
    ///
    /// For [`ContigScheme::Bare`] and [`ContigScheme::Ucsc`], delegates to
    /// [`from_name`](Self::from_name) (which handles both).  For
    /// [`ContigScheme::RefSeq`] and [`ContigScheme::GenBank`], uses
    /// assembly-specific accession lookup with strict version matching.
    ///
    /// In all cases, tries an exact match first, then a case-insensitive
    /// fallback (since VCFs in the wild use `Chr1`, `CHR1`, `chr1`, etc.).
    pub fn resolve(name: &str, scheme: ContigScheme, assembly: Assembly) -> Option<Self> {
        match scheme {
            ContigScheme::Bare | ContigScheme::Ucsc => {
                Self::from_name(name).or_else(|| Self::from_name(&name.to_ascii_lowercase()))
            }
            ContigScheme::RefSeq => Self::from_refseq(name, assembly)
                .or_else(|| Self::from_refseq(&name.to_ascii_uppercase(), assembly)),
            ContigScheme::GenBank => Self::from_genbank(name, assembly)
                .or_else(|| Self::from_genbank(&name.to_ascii_uppercase(), assembly)),
        }
    }

    /// Convert from the 5-bit field value.
    pub fn from_bits(bits: u8) -> Option<Self> {
        if bits <= Self::MAX {
            Some(ChrIndex(bits))
        } else {
            None
        }
    }

    /// Convert to the 5-bit field value.
    pub fn to_bits(self) -> u8 {
        self.0
    }
}

impl fmt::Display for ChrIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.to_name() {
            Some(name) => write!(f, "{}", name),
            None => write!(f, "?{}", self.0),
        }
    }
}

/// Canonical chromosome names indexed by ChrIndex value.
pub const CHROMOSOME_NAMES: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "M",
];

/// UCSC-style chromosome names indexed by ChrIndex value.
///
/// These match the contig names used in UCSC `.2bit` reference genome files
/// (e.g. `hg38.2bit`, `hg19.2bit`).  Note: UCSC uses `"chrM"` (not `"chrMT"`).
pub const UCSC_CHROMOSOME_NAMES: &[&str] = &[
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY", "chrM",
];

/// GRCh38 chromosome lengths and accessions.
pub const GRCH38_CHROMOSOMES: &[ChromosomeInfo] = &[
    ChromosomeInfo {
        name: "1",
        length: 248956422,
        genbank: Some("CM000663.2"),
        refseq: "NC_000001.11",
    },
    ChromosomeInfo {
        name: "2",
        length: 242193529,
        genbank: Some("CM000664.2"),
        refseq: "NC_000002.12",
    },
    ChromosomeInfo {
        name: "3",
        length: 198295559,
        genbank: Some("CM000665.2"),
        refseq: "NC_000003.12",
    },
    ChromosomeInfo {
        name: "4",
        length: 190214555,
        genbank: Some("CM000666.2"),
        refseq: "NC_000004.12",
    },
    ChromosomeInfo {
        name: "5",
        length: 181538259,
        genbank: Some("CM000667.2"),
        refseq: "NC_000005.10",
    },
    ChromosomeInfo {
        name: "6",
        length: 170805979,
        genbank: Some("CM000668.2"),
        refseq: "NC_000006.12",
    },
    ChromosomeInfo {
        name: "7",
        length: 159345973,
        genbank: Some("CM000669.2"),
        refseq: "NC_000007.14",
    },
    ChromosomeInfo {
        name: "8",
        length: 145138636,
        genbank: Some("CM000670.2"),
        refseq: "NC_000008.11",
    },
    ChromosomeInfo {
        name: "9",
        length: 138394717,
        genbank: Some("CM000671.2"),
        refseq: "NC_000009.12",
    },
    ChromosomeInfo {
        name: "10",
        length: 133797422,
        genbank: Some("CM000672.2"),
        refseq: "NC_000010.11",
    },
    ChromosomeInfo {
        name: "11",
        length: 135086622,
        genbank: Some("CM000673.2"),
        refseq: "NC_000011.10",
    },
    ChromosomeInfo {
        name: "12",
        length: 133275309,
        genbank: Some("CM000674.2"),
        refseq: "NC_000012.12",
    },
    ChromosomeInfo {
        name: "13",
        length: 114364328,
        genbank: Some("CM000675.2"),
        refseq: "NC_000013.11",
    },
    ChromosomeInfo {
        name: "14",
        length: 107043718,
        genbank: Some("CM000676.2"),
        refseq: "NC_000014.9",
    },
    ChromosomeInfo {
        name: "15",
        length: 101991189,
        genbank: Some("CM000677.2"),
        refseq: "NC_000015.10",
    },
    ChromosomeInfo {
        name: "16",
        length: 90338345,
        genbank: Some("CM000678.2"),
        refseq: "NC_000016.10",
    },
    ChromosomeInfo {
        name: "17",
        length: 83257441,
        genbank: Some("CM000679.2"),
        refseq: "NC_000017.11",
    },
    ChromosomeInfo {
        name: "18",
        length: 80373285,
        genbank: Some("CM000680.2"),
        refseq: "NC_000018.10",
    },
    ChromosomeInfo {
        name: "19",
        length: 58617616,
        genbank: Some("CM000681.2"),
        refseq: "NC_000019.10",
    },
    ChromosomeInfo {
        name: "20",
        length: 64444167,
        genbank: Some("CM000682.2"),
        refseq: "NC_000020.11",
    },
    ChromosomeInfo {
        name: "21",
        length: 46709983,
        genbank: Some("CM000683.2"),
        refseq: "NC_000021.9",
    },
    ChromosomeInfo {
        name: "22",
        length: 50818468,
        genbank: Some("CM000684.2"),
        refseq: "NC_000022.11",
    },
    ChromosomeInfo {
        name: "X",
        length: 156040895,
        genbank: Some("CM000685.2"),
        refseq: "NC_000023.11",
    },
    ChromosomeInfo {
        name: "Y",
        length: 57227415,
        genbank: Some("CM000686.2"),
        refseq: "NC_000024.10",
    },
    ChromosomeInfo {
        name: "M",
        length: 16569,
        genbank: None,
        refseq: "NC_012920.1",
    },
];

/// GRCh37 chromosome lengths and accessions.
pub const GRCH37_CHROMOSOMES: &[ChromosomeInfo] = &[
    ChromosomeInfo {
        name: "1",
        length: 249250621,
        genbank: Some("CM000663.1"),
        refseq: "NC_000001.10",
    },
    ChromosomeInfo {
        name: "2",
        length: 243199373,
        genbank: Some("CM000664.1"),
        refseq: "NC_000002.11",
    },
    ChromosomeInfo {
        name: "3",
        length: 198022430,
        genbank: Some("CM000665.1"),
        refseq: "NC_000003.11",
    },
    ChromosomeInfo {
        name: "4",
        length: 191154276,
        genbank: Some("CM000666.1"),
        refseq: "NC_000004.11",
    },
    ChromosomeInfo {
        name: "5",
        length: 180915260,
        genbank: Some("CM000667.1"),
        refseq: "NC_000005.9",
    },
    ChromosomeInfo {
        name: "6",
        length: 171115067,
        genbank: Some("CM000668.1"),
        refseq: "NC_000006.11",
    },
    ChromosomeInfo {
        name: "7",
        length: 159138663,
        genbank: Some("CM000669.1"),
        refseq: "NC_000007.13",
    },
    ChromosomeInfo {
        name: "8",
        length: 146364022,
        genbank: Some("CM000670.1"),
        refseq: "NC_000008.10",
    },
    ChromosomeInfo {
        name: "9",
        length: 141213431,
        genbank: Some("CM000671.1"),
        refseq: "NC_000009.11",
    },
    ChromosomeInfo {
        name: "10",
        length: 135534747,
        genbank: Some("CM000672.1"),
        refseq: "NC_000010.10",
    },
    ChromosomeInfo {
        name: "11",
        length: 135006516,
        genbank: Some("CM000673.1"),
        refseq: "NC_000011.9",
    },
    ChromosomeInfo {
        name: "12",
        length: 133851895,
        genbank: Some("CM000674.1"),
        refseq: "NC_000012.11",
    },
    ChromosomeInfo {
        name: "13",
        length: 115169878,
        genbank: Some("CM000675.1"),
        refseq: "NC_000013.10",
    },
    ChromosomeInfo {
        name: "14",
        length: 107349540,
        genbank: Some("CM000676.1"),
        refseq: "NC_000014.8",
    },
    ChromosomeInfo {
        name: "15",
        length: 102531392,
        genbank: Some("CM000677.1"),
        refseq: "NC_000015.9",
    },
    ChromosomeInfo {
        name: "16",
        length: 90354753,
        genbank: Some("CM000678.1"),
        refseq: "NC_000016.9",
    },
    ChromosomeInfo {
        name: "17",
        length: 81195210,
        genbank: Some("CM000679.1"),
        refseq: "NC_000017.10",
    },
    ChromosomeInfo {
        name: "18",
        length: 78077248,
        genbank: Some("CM000680.1"),
        refseq: "NC_000018.9",
    },
    ChromosomeInfo {
        name: "19",
        length: 59128983,
        genbank: Some("CM000681.1"),
        refseq: "NC_000019.9",
    },
    ChromosomeInfo {
        name: "20",
        length: 63025520,
        genbank: Some("CM000682.1"),
        refseq: "NC_000020.10",
    },
    ChromosomeInfo {
        name: "21",
        length: 48129895,
        genbank: Some("CM000683.1"),
        refseq: "NC_000021.8",
    },
    ChromosomeInfo {
        name: "22",
        length: 51304566,
        genbank: Some("CM000684.1"),
        refseq: "NC_000022.10",
    },
    ChromosomeInfo {
        name: "X",
        length: 155270560,
        genbank: Some("CM000685.1"),
        refseq: "NC_000023.10",
    },
    ChromosomeInfo {
        name: "Y",
        length: 59373566,
        genbank: Some("CM000686.1"),
        refseq: "NC_000024.9",
    },
    ChromosomeInfo {
        name: "M",
        length: 16569,
        genbank: None,
        refseq: "NC_012920.1",
    },
];

/// Get the chromosome info array for a given assembly.
pub fn chromosomes(assembly: Assembly) -> &'static [ChromosomeInfo] {
    match assembly {
        Assembly::GRCh37 => GRCH37_CHROMOSOMES,
        Assembly::GRCh38 => GRCH38_CHROMOSOMES,
    }
}

/// Validate that a position is within range for the given chromosome and assembly.
pub fn validate_position(chr: ChrIndex, pos: u32, assembly: Assembly) -> bool {
    let chroms = chromosomes(assembly);
    if let Some(info) = chroms.get(chr.0 as usize) {
        pos >= 1 && pos <= info.length
    } else {
        false
    }
}

/// GRCh38 linearized chromosome offsets.
///
/// The offset for chromosome i is the cumulative sum of all chromosome lengths
/// before it. The linearized position = `GRCH38_OFFSETS[chr_index] + pos`.
pub const GRCH38_OFFSETS: [u32; 25] = [
    0,          // chr1
    248956422,  // chr2
    491149951,  // chr3
    689445510,  // chr4
    879660065,  // chr5
    1061198324, // chr6
    1232004303, // chr7
    1391350276, // chr8
    1536488912, // chr9
    1674883629, // chr10
    1808681051, // chr11
    1943767673, // chr12
    2077042982, // chr13
    2191407310, // chr14
    2298451028, // chr15
    2400442217, // chr16
    2490780562, // chr17
    2574038003, // chr18
    2654411288, // chr19
    2713028904, // chr20
    2777473071, // chr21
    2824183054, // chr22
    2875001522, // chrX
    3031042417, // chrY
    3088269832, // chrM
];

/// GRCh37 linearized chromosome offsets.
pub const GRCH37_OFFSETS: [u32; 25] = [
    0,          // chr1
    249250621,  // chr2
    492449994,  // chr3
    690472424,  // chr4
    881626700,  // chr5
    1062541960, // chr6
    1233657027, // chr7
    1392795690, // chr8
    1539159712, // chr9
    1680373143, // chr10
    1815907890, // chr11
    1950914406, // chr12
    2084766301, // chr13
    2199936179, // chr14
    2307285719, // chr15
    2409817111, // chr16
    2500171864, // chr17
    2581367074, // chr18
    2659444322, // chr19
    2718573305, // chr20
    2781598825, // chr21
    2829728720, // chr22
    2881033286, // chrX
    3036303846, // chrY
    3095677412, // chrM
];

/// Get the linearized offset table for a given assembly.
pub fn offsets(assembly: Assembly) -> &'static [u32; 25] {
    match assembly {
        Assembly::GRCh37 => &GRCH37_OFFSETS,
        Assembly::GRCh38 => &GRCH38_OFFSETS,
    }
}

/// Convert a chromosome name and 1-based position to a linearized genome coordinate.
///
/// The linearized coordinate concatenates all chromosomes into a single
/// coordinate space, matching the encoding used in the original UVID 64-bit design.
///
/// Returns `None` if the chromosome is invalid or the position is out of range.
///
/// ```
/// use uvid::assembly::{Assembly, ChrIndex, chr_pos_to_linear};
/// assert_eq!(chr_pos_to_linear(ChrIndex(0), 1, Assembly::GRCh38), Some(1));
/// assert_eq!(chr_pos_to_linear(ChrIndex(2), 1, Assembly::GRCh38), Some(491149952));
/// ```
pub fn chr_pos_to_linear(chr: ChrIndex, pos: u32, assembly: Assembly) -> Option<u32> {
    if chr.0 > ChrIndex::MAX {
        return None;
    }
    if pos == 0 {
        return None;
    }
    let chroms = chromosomes(assembly);
    let info = chroms.get(chr.0 as usize)?;
    if pos > info.length {
        return None;
    }
    let table = offsets(assembly);
    // Checked add to ensure no overflow (though it fits in u32 for both assemblies)
    table[chr.0 as usize].checked_add(pos)
}

/// Convert a linearized genome coordinate back to a chromosome and 1-based position.
///
/// Returns `None` if the offset is 0 or beyond the end of the genome.
///
/// ```
/// use uvid::assembly::{Assembly, ChrIndex, linear_to_chr_pos};
/// assert_eq!(linear_to_chr_pos(1, Assembly::GRCh38), Some((ChrIndex(0), 1)));
/// assert_eq!(linear_to_chr_pos(491149952, Assembly::GRCh38), Some((ChrIndex(2), 1)));
/// ```
pub fn linear_to_chr_pos(linear: u32, assembly: Assembly) -> Option<(ChrIndex, u32)> {
    if linear == 0 {
        return None;
    }
    let table = offsets(assembly);
    let chroms = chromosomes(assembly);

    // Binary search: find the last chromosome whose offset is < linear
    // partition_point returns the first index where offset >= linear
    let idx = table.partition_point(|&offset| offset < linear);
    // idx is the first chromosome whose offset >= linear, so we want idx - 1
    if idx == 0 {
        return None;
    }
    let chr_idx = idx - 1;
    let pos = linear - table[chr_idx];

    // Validate the position is within chromosome bounds
    if chr_idx < chroms.len() && pos >= 1 && pos <= chroms[chr_idx].length {
        Some((ChrIndex(chr_idx as u8), pos))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chr_index_from_name() {
        assert_eq!(ChrIndex::from_name("1"), Some(ChrIndex(0)));
        assert_eq!(ChrIndex::from_name("22"), Some(ChrIndex(21)));
        assert_eq!(ChrIndex::from_name("X"), Some(ChrIndex(22)));
        assert_eq!(ChrIndex::from_name("Y"), Some(ChrIndex(23)));
        assert_eq!(ChrIndex::from_name("M"), Some(ChrIndex(24)));
        assert_eq!(ChrIndex::from_name("MT"), Some(ChrIndex(24)));
        // With chr prefix
        assert_eq!(ChrIndex::from_name("chr1"), Some(ChrIndex(0)));
        assert_eq!(ChrIndex::from_name("chrX"), Some(ChrIndex(22)));
        assert_eq!(ChrIndex::from_name("chrM"), Some(ChrIndex(24)));
        assert_eq!(ChrIndex::from_name("chrMT"), Some(ChrIndex(24)));
        // Invalid
        assert_eq!(ChrIndex::from_name("0"), None);
        assert_eq!(ChrIndex::from_name("23"), None);
        assert_eq!(ChrIndex::from_name("W"), None);
    }

    #[test]
    fn test_chr_index_roundtrip() {
        for i in 0..=ChrIndex::MAX {
            let idx = ChrIndex(i);
            let name = idx.to_name().unwrap();
            let recovered = ChrIndex::from_name(name).unwrap();
            assert_eq!(recovered, idx);
        }
    }

    #[test]
    fn test_assembly_from_str() {
        assert_eq!("GRCh37".parse::<Assembly>(), Ok(Assembly::GRCh37));
        assert_eq!("hg19".parse::<Assembly>(), Ok(Assembly::GRCh37));
        assert_eq!("37".parse::<Assembly>(), Ok(Assembly::GRCh37));
        assert_eq!("GRCh38".parse::<Assembly>(), Ok(Assembly::GRCh38));
        assert_eq!("hg38".parse::<Assembly>(), Ok(Assembly::GRCh38));
        assert_eq!("38".parse::<Assembly>(), Ok(Assembly::GRCh38));
    }

    #[test]
    fn test_assembly_bits_roundtrip() {
        for asm in [Assembly::GRCh37, Assembly::GRCh38] {
            let bits = asm.to_bits();
            let recovered = Assembly::from_bits(bits).unwrap();
            assert_eq!(recovered, asm);
        }
    }

    #[test]
    fn test_validate_position() {
        // chr1 GRCh38 length is 248956422
        assert!(validate_position(ChrIndex(0), 1, Assembly::GRCh38));
        assert!(validate_position(ChrIndex(0), 248956422, Assembly::GRCh38));
        assert!(!validate_position(ChrIndex(0), 0, Assembly::GRCh38));
        assert!(!validate_position(ChrIndex(0), 248956423, Assembly::GRCh38));
        // chrM is short
        assert!(validate_position(ChrIndex(24), 16569, Assembly::GRCh38));
        assert!(!validate_position(ChrIndex(24), 16570, Assembly::GRCh38));
    }

    #[test]
    fn test_chromosome_data_consistency() {
        // Both assemblies should have 25 chromosomes
        assert_eq!(GRCH38_CHROMOSOMES.len(), 25);
        assert_eq!(GRCH37_CHROMOSOMES.len(), 25);
        // Mitochondrial genome is the same in both
        assert_eq!(GRCH38_CHROMOSOMES[24].length, 16569);
        assert_eq!(GRCH37_CHROMOSOMES[24].length, 16569);
    }

    #[test]
    fn test_offsets_consistency() {
        // Verify offset tables are consistent with chromosome lengths
        for assembly in [Assembly::GRCh37, Assembly::GRCh38] {
            let table = offsets(assembly);
            let chroms = chromosomes(assembly);
            // chr1 always starts at offset 0
            assert_eq!(table[0], 0);
            // Each subsequent offset = previous offset + previous chromosome length
            for i in 1..25 {
                assert_eq!(
                    table[i],
                    table[i - 1] + chroms[i - 1].length,
                    "Offset mismatch at chromosome index {} for {:?}",
                    i,
                    assembly
                );
            }
        }
    }

    #[test]
    fn test_chr_pos_to_linear_chr1_start() {
        // chr1:1 in GRCh38 = offset 0 + 1 = 1
        assert_eq!(chr_pos_to_linear(ChrIndex(0), 1, Assembly::GRCh38), Some(1));
        // chr1:1 in GRCh37 = offset 0 + 1 = 1
        assert_eq!(chr_pos_to_linear(ChrIndex(0), 1, Assembly::GRCh37), Some(1));
    }

    #[test]
    fn test_chr_pos_to_linear_chr3_start() {
        // chr3:1 in GRCh38 = 491149951 + 1 = 491149952
        assert_eq!(
            chr_pos_to_linear(ChrIndex(2), 1, Assembly::GRCh38),
            Some(491149952)
        );
        // chr3:1 in GRCh37 = 492449994 + 1 = 492449995
        assert_eq!(
            chr_pos_to_linear(ChrIndex(2), 1, Assembly::GRCh37),
            Some(492449995)
        );
    }

    #[test]
    fn test_chr_pos_to_linear_chrm_end() {
        // chrM last base in GRCh38 = 3088269832 + 16569 = 3088286401
        assert_eq!(
            chr_pos_to_linear(ChrIndex(24), 16569, Assembly::GRCh38),
            Some(3088286401)
        );
    }

    #[test]
    fn test_chr_pos_to_linear_invalid() {
        // pos=0 is invalid
        assert_eq!(chr_pos_to_linear(ChrIndex(0), 0, Assembly::GRCh38), None);
        // pos beyond chromosome length
        assert_eq!(
            chr_pos_to_linear(ChrIndex(0), 248956423, Assembly::GRCh38),
            None
        );
        // Invalid chromosome index
        assert_eq!(chr_pos_to_linear(ChrIndex(25), 1, Assembly::GRCh38), None);
    }

    #[test]
    fn test_linear_to_chr_pos_chr1_start() {
        assert_eq!(
            linear_to_chr_pos(1, Assembly::GRCh38),
            Some((ChrIndex(0), 1))
        );
    }

    #[test]
    fn test_linear_to_chr_pos_chr3_start() {
        assert_eq!(
            linear_to_chr_pos(491149952, Assembly::GRCh38),
            Some((ChrIndex(2), 1))
        );
    }

    #[test]
    fn test_linear_to_chr_pos_chrm_end() {
        assert_eq!(
            linear_to_chr_pos(3088286401, Assembly::GRCh38),
            Some((ChrIndex(24), 16569))
        );
    }

    #[test]
    fn test_linear_to_chr_pos_invalid() {
        // 0 is not a valid linearized position
        assert_eq!(linear_to_chr_pos(0, Assembly::GRCh38), None);
        // Beyond end of genome
        assert_eq!(linear_to_chr_pos(3088286402, Assembly::GRCh38), None);
    }

    #[test]
    fn test_linear_roundtrip_all_chromosomes() {
        // Roundtrip for the first base of every chromosome in both assemblies
        for assembly in [Assembly::GRCh37, Assembly::GRCh38] {
            for chr_idx in 0..=ChrIndex::MAX {
                let chr = ChrIndex(chr_idx);
                let linear = chr_pos_to_linear(chr, 1, assembly).unwrap_or_else(|| {
                    panic!("Failed to linearize chr{}:1 in {:?}", chr_idx, assembly)
                });
                let (recovered_chr, recovered_pos) = linear_to_chr_pos(linear, assembly)
                    .unwrap_or_else(|| {
                        panic!("Failed to delinearize {} in {:?}", linear, assembly)
                    });
                assert_eq!(
                    recovered_chr, chr,
                    "chr mismatch for index {} in {:?}",
                    chr_idx, assembly
                );
                assert_eq!(
                    recovered_pos, 1,
                    "pos mismatch for index {} in {:?}",
                    chr_idx, assembly
                );
            }
        }
    }

    #[test]
    fn test_linear_roundtrip_last_bases() {
        // Roundtrip for the last base of every chromosome
        for assembly in [Assembly::GRCh37, Assembly::GRCh38] {
            let chroms = chromosomes(assembly);
            for chr_idx in 0..=ChrIndex::MAX {
                let chr = ChrIndex(chr_idx);
                let last_pos = chroms[chr_idx as usize].length;
                let linear = chr_pos_to_linear(chr, last_pos, assembly).unwrap();
                let (recovered_chr, recovered_pos) = linear_to_chr_pos(linear, assembly).unwrap();
                assert_eq!(recovered_chr, chr);
                assert_eq!(recovered_pos, last_pos);
            }
        }
    }

    #[test]
    fn test_linear_chromosome_boundaries() {
        // The last base of chr1 and first base of chr2 should be adjacent
        let chr1_end = chr_pos_to_linear(ChrIndex(0), 248956422, Assembly::GRCh38).unwrap();
        let chr2_start = chr_pos_to_linear(ChrIndex(1), 1, Assembly::GRCh38).unwrap();
        assert_eq!(chr2_start, chr1_end + 1);
    }

    #[test]
    fn test_linear_different_assemblies_differ() {
        // Same chromosome and position should give different linear coordinates
        // for different assemblies (because chromosome lengths differ)
        let grch38 = chr_pos_to_linear(ChrIndex(2), 1, Assembly::GRCh38).unwrap();
        let grch37 = chr_pos_to_linear(ChrIndex(2), 1, Assembly::GRCh37).unwrap();
        assert_ne!(grch38, grch37);
        // But chr1:1 is always 1 for both
        assert_eq!(
            chr_pos_to_linear(ChrIndex(0), 1, Assembly::GRCh38).unwrap(),
            1
        );
        assert_eq!(
            chr_pos_to_linear(ChrIndex(0), 1, Assembly::GRCh37).unwrap(),
            1
        );
    }

    #[test]
    fn test_linear_fits_u32() {
        // Maximum linearized position must fit in u32 for both assemblies
        for assembly in [Assembly::GRCh37, Assembly::GRCh38] {
            let chroms = chromosomes(assembly);
            let max_linear = chr_pos_to_linear(ChrIndex(24), chroms[24].length, assembly).unwrap();
            assert!(
                max_linear <= u32::MAX,
                "Max linear position exceeds u32 for {:?}",
                assembly
            );
        }
    }

    // -------------------------------------------------------------------
    // ContigScheme, classify_chrom, detect_contig_scheme
    // -------------------------------------------------------------------

    #[test]
    fn test_classify_chrom_bare() {
        assert_eq!(classify_chrom("1"), ContigScheme::Bare);
        assert_eq!(classify_chrom("22"), ContigScheme::Bare);
        assert_eq!(classify_chrom("X"), ContigScheme::Bare);
        assert_eq!(classify_chrom("Y"), ContigScheme::Bare);
        assert_eq!(classify_chrom("MT"), ContigScheme::Bare);
    }

    #[test]
    fn test_classify_chrom_ucsc() {
        assert_eq!(classify_chrom("chr1"), ContigScheme::Ucsc);
        assert_eq!(classify_chrom("chrX"), ContigScheme::Ucsc);
        assert_eq!(classify_chrom("chrM"), ContigScheme::Ucsc);
        assert_eq!(classify_chrom("chrMT"), ContigScheme::Ucsc);
        // Case variants
        assert_eq!(classify_chrom("Chr1"), ContigScheme::Ucsc);
        assert_eq!(classify_chrom("CHR1"), ContigScheme::Ucsc);
        assert_eq!(classify_chrom("ChrX"), ContigScheme::Ucsc);
    }

    #[test]
    fn test_classify_chrom_refseq() {
        assert_eq!(classify_chrom("NC_000001.11"), ContigScheme::RefSeq);
        assert_eq!(classify_chrom("NC_000023.11"), ContigScheme::RefSeq);
        assert_eq!(classify_chrom("NC_012920.1"), ContigScheme::RefSeq);
        // Case variant (unlikely but tolerated)
        assert_eq!(classify_chrom("nc_000001.11"), ContigScheme::RefSeq);
    }

    #[test]
    fn test_classify_chrom_genbank() {
        assert_eq!(classify_chrom("CM000663.2"), ContigScheme::GenBank);
        assert_eq!(classify_chrom("CM000686.1"), ContigScheme::GenBank);
        // Case variant
        assert_eq!(classify_chrom("cm000663.2"), ContigScheme::GenBank);
    }

    #[test]
    fn test_detect_contig_scheme_empty() {
        assert_eq!(detect_contig_scheme(&[]), ContigScheme::Bare);
    }

    #[test]
    fn test_detect_contig_scheme_bare() {
        assert_eq!(
            detect_contig_scheme(&["1", "2", "3", "X", "Y", "MT"]),
            ContigScheme::Bare
        );
    }

    #[test]
    fn test_detect_contig_scheme_ucsc() {
        assert_eq!(
            detect_contig_scheme(&["chr1", "chr2", "chrX"]),
            ContigScheme::Ucsc
        );
    }

    #[test]
    fn test_detect_contig_scheme_refseq() {
        assert_eq!(
            detect_contig_scheme(&["NC_000001.11", "NC_000002.12"]),
            ContigScheme::RefSeq
        );
    }

    #[test]
    fn test_detect_contig_scheme_genbank() {
        assert_eq!(
            detect_contig_scheme(&["CM000663.2", "CM000664.2"]),
            ContigScheme::GenBank
        );
    }

    #[test]
    fn test_detect_contig_scheme_mixed_prefers_specific() {
        // If both bare and UCSC-like contigs appear, UCSC wins
        assert_eq!(
            detect_contig_scheme(&["chr1", "some_scaffold"]),
            ContigScheme::Ucsc
        );
        // If RefSeq appears alongside bare, RefSeq wins
        assert_eq!(
            detect_contig_scheme(&["1", "NC_000001.11"]),
            ContigScheme::RefSeq
        );
    }

    // -------------------------------------------------------------------
    // ChrIndex::from_refseq
    // -------------------------------------------------------------------

    #[test]
    fn test_from_refseq_grch38() {
        assert_eq!(
            ChrIndex::from_refseq("NC_000001.11", Assembly::GRCh38),
            Some(ChrIndex(0))
        );
        assert_eq!(
            ChrIndex::from_refseq("NC_000022.11", Assembly::GRCh38),
            Some(ChrIndex(21))
        );
        assert_eq!(
            ChrIndex::from_refseq("NC_000023.11", Assembly::GRCh38),
            Some(ChrIndex(22)), // chrX
        );
        assert_eq!(
            ChrIndex::from_refseq("NC_000024.10", Assembly::GRCh38),
            Some(ChrIndex(23)), // chrY
        );
        assert_eq!(
            ChrIndex::from_refseq("NC_012920.1", Assembly::GRCh38),
            Some(ChrIndex(24)), // chrM
        );
    }

    #[test]
    fn test_from_refseq_grch37() {
        assert_eq!(
            ChrIndex::from_refseq("NC_000001.10", Assembly::GRCh37),
            Some(ChrIndex(0))
        );
        assert_eq!(
            ChrIndex::from_refseq("NC_000023.10", Assembly::GRCh37),
            Some(ChrIndex(22)), // chrX
        );
        // chrM is the same accession+version in both assemblies
        assert_eq!(
            ChrIndex::from_refseq("NC_012920.1", Assembly::GRCh37),
            Some(ChrIndex(24)),
        );
    }

    #[test]
    fn test_from_refseq_version_mismatch() {
        // GRCh37 version for chr1 should NOT match GRCh38
        assert_eq!(
            ChrIndex::from_refseq("NC_000001.10", Assembly::GRCh38),
            None
        );
        // GRCh38 version for chr1 should NOT match GRCh37
        assert_eq!(
            ChrIndex::from_refseq("NC_000001.11", Assembly::GRCh37),
            None
        );
    }

    #[test]
    fn test_from_refseq_invalid() {
        assert_eq!(ChrIndex::from_refseq("NC_999999.1", Assembly::GRCh38), None);
        assert_eq!(ChrIndex::from_refseq("chr1", Assembly::GRCh38), None);
        assert_eq!(ChrIndex::from_refseq("", Assembly::GRCh38), None);
    }

    // -------------------------------------------------------------------
    // ChrIndex::from_genbank
    // -------------------------------------------------------------------

    #[test]
    fn test_from_genbank_grch38() {
        assert_eq!(
            ChrIndex::from_genbank("CM000663.2", Assembly::GRCh38),
            Some(ChrIndex(0))
        );
        assert_eq!(
            ChrIndex::from_genbank("CM000685.2", Assembly::GRCh38),
            Some(ChrIndex(22)), // chrX
        );
        assert_eq!(
            ChrIndex::from_genbank("CM000686.2", Assembly::GRCh38),
            Some(ChrIndex(23)), // chrY
        );
    }

    #[test]
    fn test_from_genbank_grch37() {
        assert_eq!(
            ChrIndex::from_genbank("CM000663.1", Assembly::GRCh37),
            Some(ChrIndex(0))
        );
    }

    #[test]
    fn test_from_genbank_version_mismatch() {
        assert_eq!(ChrIndex::from_genbank("CM000663.1", Assembly::GRCh38), None);
        assert_eq!(ChrIndex::from_genbank("CM000663.2", Assembly::GRCh37), None);
    }

    #[test]
    fn test_from_genbank_chrm_has_no_genbank() {
        // chrM has genbank: None in both assemblies
        assert_eq!(
            ChrIndex::from_genbank("NC_012920.1", Assembly::GRCh38),
            None
        );
    }

    // -------------------------------------------------------------------
    // ChrIndex::to_ucsc_name
    // -------------------------------------------------------------------

    #[test]
    fn test_to_ucsc_name() {
        assert_eq!(ChrIndex(0).to_ucsc_name(), Some("chr1"));
        assert_eq!(ChrIndex(21).to_ucsc_name(), Some("chr22"));
        assert_eq!(ChrIndex(22).to_ucsc_name(), Some("chrX"));
        assert_eq!(ChrIndex(23).to_ucsc_name(), Some("chrY"));
        assert_eq!(ChrIndex(24).to_ucsc_name(), Some("chrM"));
        assert_eq!(ChrIndex(25).to_ucsc_name(), None);
    }

    #[test]
    fn test_ucsc_name_roundtrip() {
        // Every UCSC name should resolve back to the same ChrIndex
        for i in 0..=ChrIndex::MAX {
            let idx = ChrIndex(i);
            let ucsc = idx.to_ucsc_name().unwrap();
            assert_eq!(
                ChrIndex::from_name(ucsc),
                Some(idx),
                "roundtrip failed for {}",
                ucsc
            );
        }
    }

    // -------------------------------------------------------------------
    // ChrIndex::resolve (scheme-aware, case-insensitive)
    // -------------------------------------------------------------------

    #[test]
    fn test_resolve_bare() {
        let asm = Assembly::GRCh38;
        assert_eq!(
            ChrIndex::resolve("1", ContigScheme::Bare, asm),
            Some(ChrIndex(0))
        );
        assert_eq!(
            ChrIndex::resolve("X", ContigScheme::Bare, asm),
            Some(ChrIndex(22))
        );
        assert_eq!(
            ChrIndex::resolve("MT", ContigScheme::Bare, asm),
            Some(ChrIndex(24))
        );
    }

    #[test]
    fn test_resolve_ucsc() {
        let asm = Assembly::GRCh38;
        assert_eq!(
            ChrIndex::resolve("chr1", ContigScheme::Ucsc, asm),
            Some(ChrIndex(0))
        );
        assert_eq!(
            ChrIndex::resolve("chrX", ContigScheme::Ucsc, asm),
            Some(ChrIndex(22))
        );
        assert_eq!(
            ChrIndex::resolve("chrM", ContigScheme::Ucsc, asm),
            Some(ChrIndex(24))
        );
    }

    #[test]
    fn test_resolve_case_insensitive_ucsc() {
        let asm = Assembly::GRCh38;
        assert_eq!(
            ChrIndex::resolve("Chr1", ContigScheme::Ucsc, asm),
            Some(ChrIndex(0))
        );
        assert_eq!(
            ChrIndex::resolve("CHR1", ContigScheme::Ucsc, asm),
            Some(ChrIndex(0))
        );
        assert_eq!(
            ChrIndex::resolve("ChrX", ContigScheme::Ucsc, asm),
            Some(ChrIndex(22))
        );
        assert_eq!(
            ChrIndex::resolve("CHRX", ContigScheme::Ucsc, asm),
            Some(ChrIndex(22))
        );
        assert_eq!(
            ChrIndex::resolve("ChrM", ContigScheme::Ucsc, asm),
            Some(ChrIndex(24))
        );
        assert_eq!(
            ChrIndex::resolve("CHRMT", ContigScheme::Ucsc, asm),
            Some(ChrIndex(24))
        );
    }

    #[test]
    fn test_resolve_case_insensitive_bare() {
        let asm = Assembly::GRCh38;
        // Lowercase x/y should work
        assert_eq!(
            ChrIndex::resolve("x", ContigScheme::Bare, asm),
            Some(ChrIndex(22))
        );
        assert_eq!(
            ChrIndex::resolve("y", ContigScheme::Bare, asm),
            Some(ChrIndex(23))
        );
        // Mt (mixed case)
        assert_eq!(
            ChrIndex::resolve("Mt", ContigScheme::Bare, asm),
            Some(ChrIndex(24))
        );
        assert_eq!(
            ChrIndex::resolve("mT", ContigScheme::Bare, asm),
            Some(ChrIndex(24))
        );
    }

    #[test]
    fn test_resolve_refseq() {
        assert_eq!(
            ChrIndex::resolve("NC_000001.11", ContigScheme::RefSeq, Assembly::GRCh38),
            Some(ChrIndex(0))
        );
        assert_eq!(
            ChrIndex::resolve("NC_000001.10", ContigScheme::RefSeq, Assembly::GRCh37),
            Some(ChrIndex(0))
        );
        // Version mismatch
        assert_eq!(
            ChrIndex::resolve("NC_000001.10", ContigScheme::RefSeq, Assembly::GRCh38),
            None
        );
    }

    #[test]
    fn test_resolve_genbank() {
        assert_eq!(
            ChrIndex::resolve("CM000663.2", ContigScheme::GenBank, Assembly::GRCh38),
            Some(ChrIndex(0))
        );
        assert_eq!(
            ChrIndex::resolve("CM000663.1", ContigScheme::GenBank, Assembly::GRCh37),
            Some(ChrIndex(0))
        );
        // Version mismatch
        assert_eq!(
            ChrIndex::resolve("CM000663.1", ContigScheme::GenBank, Assembly::GRCh38),
            None
        );
    }

    #[test]
    fn test_resolve_refseq_all_grch38_chromosomes() {
        let asm = Assembly::GRCh38;
        let chroms = chromosomes(asm);
        for (i, info) in chroms.iter().enumerate() {
            let result = ChrIndex::resolve(info.refseq, ContigScheme::RefSeq, asm);
            assert_eq!(
                result,
                Some(ChrIndex(i as u8)),
                "RefSeq {} should resolve to ChrIndex({})",
                info.refseq,
                i,
            );
        }
    }

    #[test]
    fn test_resolve_genbank_all_grch38_chromosomes() {
        let asm = Assembly::GRCh38;
        let chroms = chromosomes(asm);
        for (i, info) in chroms.iter().enumerate() {
            if let Some(gb) = info.genbank {
                let result = ChrIndex::resolve(gb, ContigScheme::GenBank, asm);
                assert_eq!(
                    result,
                    Some(ChrIndex(i as u8)),
                    "GenBank {} should resolve to ChrIndex({})",
                    gb,
                    i,
                );
            }
        }
    }
}
