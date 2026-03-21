/// Genome assembly definitions for GRCh37 and GRCh38.
///
/// Provides chromosome metadata and encoding/decoding between
/// chromosome name + position and the compact (chr_index, position) form
/// used in the UVID 128-bit layout.
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
}
