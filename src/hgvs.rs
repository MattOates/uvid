//! HGVS genomic notation parsing and formatting.
//!
//! Supports the `g.` (genomic) and `m.` (mitochondrial) coordinate systems
//! from the [HGVS Nomenclature](https://hgvs-nomenclature.org/).  Other
//! coordinate systems (`c.`, `n.`, `p.`, `r.`) are intentionally unsupported;
//! use the [`ferro-hgvs`](https://crates.io/crates/hgvs) crate if you need
//! transcript-level notation.
//!
//! # HGVS ↔ VCF coordinate differences
//!
//! HGVS and VCF represent the same variants differently:
//!
//! | Type         | HGVS example                    | VCF (POS, REF, ALT)       |
//! |--------------|---------------------------------|---------------------------|
//! | Substitution | `NC_000001.11:g.12345A>G`       | 12345, A, G               |
//! | Insertion    | `NC_000001.11:g.12345_12346insACT` | 12345, \<anchor\>, \<anchor\>ACT |
//! | Deletion     | `NC_000001.11:g.12345_12347del` | 12344, \<anchor\>XXX, \<anchor\> |
//! | Delins       | `NC_000001.11:g.12345_12347delinsACT` | 12345, XXX, ACT     |
//! | Duplication  | `NC_000001.11:g.12345_12347dup` | (insertion of dup seq)    |
//! | Inversion    | `NC_000001.11:g.12345_12347inv` | 12345, XXX, rev-comp(XXX) |
//!
//! Insertions and deletions require a reference genome to supply the VCF
//! "anchor base" and deleted sequence.  Only pure substitutions are
//! reference-free.
//!
//! # Assembly auto-detection
//!
//! The assembly is inferred from the RefSeq accession version number:
//! `NC_000001.11` → GRCh38, `NC_000001.10` → GRCh37.  An optional
//! `expected_assembly` parameter validates the inferred assembly.

use std::fmt;

use crate::assembly::{self, Assembly, ChrIndex};
use crate::normalize::ReferenceGenome;
use crate::uvid128::{Uvid128, Variant};

// ── Error types ────────────────────────────────────────────────────────────

/// Errors that can occur during HGVS parsing or conversion.
#[derive(Debug)]
pub enum HgvsError {
    /// The HGVS string could not be parsed.
    Parse(String),

    /// The coordinate system is not supported (c., n., p., r.).
    /// We only support g. (genomic) and m. (mitochondrial).
    UnsupportedCoordinateSystem { system: char, hint: &'static str },

    /// The RefSeq accession could not be resolved to a chromosome.
    UnknownAccession(String),

    /// The assembly inferred from the accession does not match the expected one.
    AssemblyMismatch {
        inferred: Assembly,
        expected: Assembly,
    },

    /// A reference genome is required but was not provided.
    ReferenceRequired { reason: &'static str },

    /// An error occurred while accessing the reference genome.
    Reference(crate::normalize::ReferenceError),

    /// The UVID could not be decoded (malformed data).
    DecodeFailed,

    /// Allele is in length-mode (>20bp), so exact sequence is unavailable.
    LengthModeWarning { field: &'static str, length: usize },
}

impl fmt::Display for HgvsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HgvsError::Parse(msg) => write!(f, "HGVS parse error: {}", msg),
            HgvsError::UnsupportedCoordinateSystem { system, hint } => {
                write!(
                    f,
                    "coordinate system '{}.'' is not supported; {}",
                    system, hint
                )
            }
            HgvsError::UnknownAccession(acc) => {
                write!(f, "unknown RefSeq accession: {}", acc)
            }
            HgvsError::AssemblyMismatch { inferred, expected } => {
                write!(
                    f,
                    "assembly mismatch: accession implies {} but expected {}",
                    inferred, expected
                )
            }
            HgvsError::ReferenceRequired { reason } => {
                write!(f, "reference genome required: {}", reason)
            }
            HgvsError::Reference(e) => write!(f, "reference error: {}", e),
            HgvsError::DecodeFailed => write!(f, "failed to decode UVID"),
            HgvsError::LengthModeWarning { field, length } => {
                write!(
                    f,
                    "{} allele ({} bp) is in length-mode; exact sequence unavailable, \
                     using N-repeats as approximation",
                    field, length
                )
            }
        }
    }
}

impl std::error::Error for HgvsError {}

impl From<crate::normalize::ReferenceError> for HgvsError {
    fn from(e: crate::normalize::ReferenceError) -> Self {
        HgvsError::Reference(e)
    }
}

// ── Types ──────────────────────────────────────────────────────────────────

/// HGVS coordinate system prefix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateSystem {
    /// Genomic (`g.`) — used with chromosomal RefSeq accessions (NC_*).
    Genomic,
    /// Mitochondrial (`m.`) — used with the mitochondrial RefSeq accession.
    Mitochondrial,
}

impl CoordinateSystem {
    /// The single-character prefix for this coordinate system.
    pub fn prefix(self) -> char {
        match self {
            CoordinateSystem::Genomic => 'g',
            CoordinateSystem::Mitochondrial => 'm',
        }
    }
}

impl fmt::Display for CoordinateSystem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}.", self.prefix())
    }
}

/// The nucleotide-level edit described by an HGVS expression.
///
/// Design borrowed from ferro-hgvs's `NaEdit` enum, simplified for
/// genomic-only notation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NaEdit {
    /// Substitution: single base change (e.g., `A>G`).
    Substitution { reference: u8, alternative: u8 },

    /// Deletion (e.g., `del`, `delACG`).
    /// `sequence` is the deleted bases if specified in the notation.
    Deletion { sequence: Option<Vec<u8>> },

    /// Insertion (e.g., `insACGT`).
    Insertion { sequence: Vec<u8> },

    /// Deletion-insertion (e.g., `delinsACGT`, `delACGinsT`).
    Delins {
        deleted: Option<Vec<u8>>,
        inserted: Vec<u8>,
    },

    /// Duplication (e.g., `dup`, `dupACG`).
    Duplication { sequence: Option<Vec<u8>> },

    /// Inversion (e.g., `inv`).
    Inversion,

    /// Identity / no change (e.g., `=`).
    Identity,
}

/// A parsed HGVS genomic variant.
///
/// ```text
/// NC_000001.11:g.12345A>G
/// ^^^^^^^^^^^^  ^ ^^^^^^^^
/// accession     |  edit
///            coord_system
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HgvsVariant {
    /// The RefSeq accession including version (e.g., `NC_000001.11`).
    pub accession: String,
    /// The coordinate system (`g.` or `m.`).
    pub coordinate_system: CoordinateSystem,
    /// 1-based start position.
    pub start: u64,
    /// 1-based end position (same as start for substitutions and single-base edits).
    pub end: u64,
    /// The nucleotide edit.
    pub edit: NaEdit,
}

impl fmt::Display for HgvsVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.accession, self.coordinate_system)?;
        match &self.edit {
            NaEdit::Substitution {
                reference,
                alternative,
            } => {
                write!(
                    f,
                    "{}{}>{}",
                    self.start, *reference as char, *alternative as char
                )
            }
            NaEdit::Deletion { sequence } => {
                if self.start == self.end {
                    write!(f, "{}del", self.start)?;
                } else {
                    write!(f, "{}_{}", self.start, self.end)?;
                    write!(f, "del")?;
                }
                if let Some(seq) = sequence {
                    write!(f, "{}", std::str::from_utf8(seq).unwrap_or("?"))?;
                }
                Ok(())
            }
            NaEdit::Insertion { sequence } => {
                write!(
                    f,
                    "{}_{}ins{}",
                    self.start,
                    self.end,
                    std::str::from_utf8(sequence).unwrap_or("?")
                )
            }
            NaEdit::Delins { deleted, inserted } => {
                if self.start == self.end {
                    write!(f, "{}", self.start)?;
                } else {
                    write!(f, "{}_{}", self.start, self.end)?;
                }
                write!(f, "del")?;
                if let Some(seq) = deleted {
                    write!(f, "{}", std::str::from_utf8(seq).unwrap_or("?"))?;
                }
                write!(f, "ins{}", std::str::from_utf8(inserted).unwrap_or("?"))
            }
            NaEdit::Duplication { sequence } => {
                if self.start == self.end {
                    write!(f, "{}dup", self.start)?;
                } else {
                    write!(f, "{}_{}dup", self.start, self.end)?;
                }
                if let Some(seq) = sequence {
                    write!(f, "{}", std::str::from_utf8(seq).unwrap_or("?"))?;
                }
                Ok(())
            }
            NaEdit::Inversion => {
                write!(f, "{}_{}inv", self.start, self.end)
            }
            NaEdit::Identity => {
                write!(f, "{}=", self.start)
            }
        }
    }
}

// ── Parser ─────────────────────────────────────────────────────────────────

/// Parse an HGVS genomic variant string.
///
/// Accepts the form: `<accession>:<prefix>.<location><edit>`
///
/// # Examples
///
/// ```
/// use uvid::hgvs::parse;
///
/// let v = parse("NC_000001.11:g.12345A>G").unwrap();
/// assert_eq!(v.start, 12345);
/// ```
pub fn parse(input: &str) -> Result<HgvsVariant, HgvsError> {
    let input = input.trim();

    // Split on ':'
    let colon_pos = input
        .find(':')
        .ok_or_else(|| HgvsError::Parse("expected ':' separating accession and location".into()))?;
    let accession = &input[..colon_pos];
    let rest = &input[colon_pos + 1..];

    // Validate accession looks like a RefSeq accession
    if !accession.starts_with("NC_") {
        return Err(HgvsError::Parse(format!(
            "expected RefSeq accession starting with 'NC_', got '{}'",
            accession
        )));
    }

    // Parse coordinate system prefix
    if rest.len() < 2 {
        return Err(HgvsError::Parse(
            "expected coordinate system prefix (e.g., 'g.')".into(),
        ));
    }
    let prefix_char = rest.as_bytes()[0] as char;
    if rest.as_bytes()[1] != b'.' {
        return Err(HgvsError::Parse(format!(
            "expected '.' after coordinate system prefix '{}', got '{}'",
            prefix_char,
            rest.as_bytes()[1] as char
        )));
    }
    let coord_sys = match prefix_char {
        'g' => CoordinateSystem::Genomic,
        'm' => CoordinateSystem::Mitochondrial,
        'c' | 'n' => {
            return Err(HgvsError::UnsupportedCoordinateSystem {
                system: prefix_char,
                hint: "use the ferro-hgvs crate for transcript-level (c./n.) notation",
            });
        }
        'p' => {
            return Err(HgvsError::UnsupportedCoordinateSystem {
                system: prefix_char,
                hint: "use the ferro-hgvs crate for protein-level (p.) notation",
            });
        }
        'r' => {
            return Err(HgvsError::UnsupportedCoordinateSystem {
                system: prefix_char,
                hint: "use the ferro-hgvs crate for RNA-level (r.) notation",
            });
        }
        _ => {
            return Err(HgvsError::Parse(format!(
                "unknown coordinate system prefix: '{}'",
                prefix_char
            )));
        }
    };

    let body = &rest[2..]; // everything after "g." or "m."

    // Parse location and edit from the body.
    // Strategy: scan for the first non-digit, non-underscore character after
    // parsing the position(s), then dispatch on the edit type.
    parse_location_and_edit(accession, coord_sys, body)
}

/// Parse the location + edit portion of an HGVS expression (after the `g.`/`m.` prefix).
fn parse_location_and_edit(
    accession: &str,
    coord_sys: CoordinateSystem,
    body: &str,
) -> Result<HgvsVariant, HgvsError> {
    // Parse the first position (always present)
    let (start, rest) = parse_u64(body)?;

    // Check for range (underscore separator)
    let (end, rest) = if let Some(after_underscore) = rest.strip_prefix('_') {
        let (end_pos, r) = parse_u64(after_underscore)?;
        (end_pos, r)
    } else {
        (start, rest)
    };

    // Dispatch on edit type
    if rest.is_empty() {
        return Err(HgvsError::Parse("expected edit type after position".into()));
    }

    // Identity: "="
    if rest == "=" {
        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start,
            end,
            edit: NaEdit::Identity,
        });
    }

    // Substitution: single base > single base (e.g., "A>G")
    if rest.len() == 3 && rest.as_bytes()[1] == b'>' {
        let ref_base = rest.as_bytes()[0].to_ascii_uppercase();
        let alt_base = rest.as_bytes()[2].to_ascii_uppercase();
        if !is_dna_base(ref_base) || !is_dna_base(alt_base) {
            return Err(HgvsError::Parse(format!(
                "invalid bases in substitution: '{}'",
                rest
            )));
        }
        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start,
            end: start, // substitutions are single-position
            edit: NaEdit::Substitution {
                reference: ref_base,
                alternative: alt_base,
            },
        });
    }

    // Deletion-insertion: "delXXXinsYYY" or "delinsYYY"
    if let Some(delins_rest) = rest.strip_prefix("del") {
        if let Some(ins_idx) = delins_rest.find("ins") {
            let del_seq_str = &delins_rest[..ins_idx];
            let ins_seq_str = &delins_rest[ins_idx + 3..];
            let deleted = if del_seq_str.is_empty() {
                None
            } else {
                Some(parse_dna_sequence(del_seq_str)?)
            };
            let inserted = parse_dna_sequence(ins_seq_str)?;
            if inserted.is_empty() {
                return Err(HgvsError::Parse(
                    "delins must have a non-empty inserted sequence".into(),
                ));
            }
            return Ok(HgvsVariant {
                accession: accession.to_string(),
                coordinate_system: coord_sys,
                start,
                end,
                edit: NaEdit::Delins { deleted, inserted },
            });
        }

        // Pure deletion: "del" or "delACG"
        let del_seq = if delins_rest.is_empty() {
            None
        } else {
            Some(parse_dna_sequence(delins_rest)?)
        };
        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start,
            end,
            edit: NaEdit::Deletion { sequence: del_seq },
        });
    }

    // Insertion: "insACGT"
    if let Some(ins_rest) = rest.strip_prefix("ins") {
        let sequence = parse_dna_sequence(ins_rest)?;
        if sequence.is_empty() {
            return Err(HgvsError::Parse(
                "insertion must have a non-empty sequence".into(),
            ));
        }
        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start,
            end,
            edit: NaEdit::Insertion { sequence },
        });
    }

    // Duplication: "dup" or "dupACG"
    if let Some(dup_rest) = rest.strip_prefix("dup") {
        let sequence = if dup_rest.is_empty() {
            None
        } else {
            Some(parse_dna_sequence(dup_rest)?)
        };
        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start,
            end,
            edit: NaEdit::Duplication { sequence },
        });
    }

    // Inversion: "inv"
    if rest == "inv" {
        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start,
            end,
            edit: NaEdit::Inversion,
        });
    }

    Err(HgvsError::Parse(format!(
        "unrecognized edit type: '{}'",
        rest
    )))
}

/// Parse a u64 from the start of a string, returning the value and the remaining string.
fn parse_u64(s: &str) -> Result<(u64, &str), HgvsError> {
    let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
    if end == 0 {
        return Err(HgvsError::Parse(format!(
            "expected a number, got '{}'",
            s.chars().next().unwrap_or('?')
        )));
    }
    let value: u64 = s[..end]
        .parse()
        .map_err(|e| HgvsError::Parse(format!("invalid number: {}", e)))?;
    Ok((value, &s[end..]))
}

/// Parse a DNA sequence string (uppercase ACGT only).
fn parse_dna_sequence(s: &str) -> Result<Vec<u8>, HgvsError> {
    let upper = s.to_ascii_uppercase();
    let bytes = upper.as_bytes();
    for &b in bytes {
        if !is_dna_base(b) {
            return Err(HgvsError::Parse(format!(
                "invalid DNA character '{}' in sequence '{}'",
                b as char, s
            )));
        }
    }
    Ok(bytes.to_vec())
}

/// Check if a byte is a valid DNA base (A, C, G, T).
#[inline]
fn is_dna_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

// ── Assembly auto-detection ────────────────────────────────────────────────

/// Detect the assembly from a RefSeq accession by trying both assemblies.
///
/// Returns the assembly and chromosome index, or an error if the accession
/// is not recognized in either assembly.
pub fn detect_assembly(accession: &str) -> Result<(Assembly, ChrIndex), HgvsError> {
    // Try GRCh38 first (more common in modern data)
    if let Some(chr) = ChrIndex::from_refseq(accession, Assembly::GRCh38) {
        return Ok((Assembly::GRCh38, chr));
    }
    if let Some(chr) = ChrIndex::from_refseq(accession, Assembly::GRCh37) {
        return Ok((Assembly::GRCh37, chr));
    }
    Err(HgvsError::UnknownAccession(accession.to_string()))
}

// ── HGVS → UVID conversion ────────────────────────────────────────────────

/// Convert an HGVS genomic variant to a UVID.
///
/// # Arguments
///
/// * `hgvs` - The HGVS string to convert.
/// * `reference` - Optional reference genome for anchor base lookups
///   (required for insertions, deletions, delins, duplications, inversions).
/// * `expected_assembly` - Optional assembly to validate against the accession.
///
/// # Returns
///
/// The UVID as a `Uvid128`, or an error if the conversion fails.
pub fn hgvs_to_uvid(
    hgvs_str: &str,
    reference: Option<&mut dyn ReferenceGenome>,
    expected_assembly: Option<Assembly>,
) -> Result<Uvid128, HgvsError> {
    let hgvs = parse(hgvs_str)?;

    // Detect assembly from accession
    let (assembly, chr) = detect_assembly(&hgvs.accession)?;

    // Validate against expected assembly if provided
    if let Some(expected) = expected_assembly {
        if assembly != expected {
            return Err(HgvsError::AssemblyMismatch {
                inferred: assembly,
                expected,
            });
        }
    }

    // Get the UCSC chromosome name for reference lookups
    let ucsc_name = chr.to_ucsc_name().ok_or(HgvsError::DecodeFailed)?;

    // Convert HGVS coordinates to VCF-style (pos, ref, alt)
    let (vcf_pos, vcf_ref, vcf_alt) = hgvs_edit_to_vcf(&hgvs, chr, assembly, ucsc_name, reference)?;

    // Encode into UVID
    Uvid128::encode(chr, vcf_pos, &vcf_ref, &vcf_alt, assembly)
        .ok_or_else(|| HgvsError::Parse("failed to encode UVID from VCF coordinates".into()))
}

/// Convert an HGVS edit to VCF-style coordinates (pos, ref_allele, alt_allele).
///
/// HGVS uses 1-based coordinates without anchor bases for indels.
/// VCF uses 1-based coordinates WITH an anchor base for indels.
fn hgvs_edit_to_vcf(
    hgvs: &HgvsVariant,
    _chr: ChrIndex,
    _assembly: Assembly,
    ucsc_name: &str,
    reference: Option<&mut dyn ReferenceGenome>,
) -> Result<(u32, Vec<u8>, Vec<u8>), HgvsError> {
    let start = hgvs.start as u32;
    let end = hgvs.end as u32;

    match &hgvs.edit {
        NaEdit::Substitution {
            reference: ref_base,
            alternative: alt_base,
        } => {
            // Direct mapping: HGVS pos = VCF pos, no anchor needed
            Ok((start, vec![*ref_base], vec![*alt_base]))
        }

        NaEdit::Deletion { sequence } => {
            // Deletion: need anchor base before the deleted region.
            // HGVS: g.start_end del  →  VCF: pos=start-1, ref=anchor+deleted, alt=anchor
            let ref_genome = reference.ok_or(HgvsError::ReferenceRequired {
                reason: "deletions need an anchor base from the reference genome",
            })?;

            let anchor_pos = start - 1; // 1-based position before deletion
                                        // Fetch anchor base (0-based coords for reference)
            let anchor = ref_genome.get_sequence(ucsc_name, anchor_pos - 1, anchor_pos)?;

            // Get deleted sequence from reference if not provided in HGVS
            let deleted = match sequence {
                Some(seq) => seq.clone(),
                None => {
                    // Fetch from reference: positions start..=end (1-based) → 0-based [start-1, end)
                    ref_genome.get_sequence(ucsc_name, start - 1, end)?
                }
            };

            let mut vcf_ref = anchor.clone();
            vcf_ref.extend_from_slice(&deleted);
            let vcf_alt = anchor;

            Ok((anchor_pos, vcf_ref, vcf_alt))
        }

        NaEdit::Insertion { sequence } => {
            // Insertion: g.start_end insXXX  (between positions start and end)
            // VCF: pos=start, ref=anchor, alt=anchor+inserted
            let ref_genome = reference.ok_or(HgvsError::ReferenceRequired {
                reason: "insertions need an anchor base from the reference genome",
            })?;

            // The anchor is the base at position 'start' (the left flanking base)
            let anchor = ref_genome.get_sequence(ucsc_name, start - 1, start)?;

            let mut vcf_alt = anchor.clone();
            vcf_alt.extend_from_slice(sequence);

            Ok((start, anchor, vcf_alt))
        }

        NaEdit::Delins { deleted, inserted } => {
            // Deletion-insertion: g.start_end delinsXXX
            // VCF: pos=start, ref=deleted_bases, alt=inserted_bases
            // (no anchor base needed — the deleted and inserted sequences
            // replace each other at the same position)
            let vcf_ref = match deleted {
                Some(seq) => seq.clone(),
                None => {
                    // Need reference to get deleted bases
                    let ref_genome = reference.ok_or(HgvsError::ReferenceRequired {
                        reason: "delins without explicit deleted sequence needs reference",
                    })?;
                    ref_genome.get_sequence(ucsc_name, start - 1, end)?
                }
            };

            Ok((start, vcf_ref, inserted.clone()))
        }

        NaEdit::Duplication { sequence } => {
            // Duplication: semantically an insertion of the duplicated sequence
            // after the duplicated region.
            // g.start_end dup → insert a copy of start..end after position end.
            // VCF: pos=end, ref=anchor, alt=anchor+dup_seq
            let ref_genome = reference.ok_or(HgvsError::ReferenceRequired {
                reason: "duplications need the reference to determine the duplicated sequence",
            })?;

            let dup_seq = match sequence {
                Some(seq) => seq.clone(),
                None => ref_genome.get_sequence(ucsc_name, start - 1, end)?,
            };

            // Anchor is the last base of the duplicated region
            let anchor = ref_genome.get_sequence(ucsc_name, end - 1, end)?;

            let mut vcf_alt = anchor.clone();
            vcf_alt.extend_from_slice(&dup_seq);

            Ok((end, anchor, vcf_alt))
        }

        NaEdit::Inversion => {
            // Inversion: ref = original, alt = reverse complement
            let ref_genome = reference.ok_or(HgvsError::ReferenceRequired {
                reason: "inversions need the reference to determine the original sequence",
            })?;

            let original = ref_genome.get_sequence(ucsc_name, start - 1, end)?;

            if original.len() > 12 {
                // Warn: long inversions will use length-mode encoding, losing sequence info.
                // This is advisory — we still encode it.
                eprintln!(
                    "warning: inversion at {}:{}.{}_{}inv is {} bp; \
                     UVID length-mode encoding (>20 bp) will lose exact sequence",
                    hgvs.accession,
                    hgvs.coordinate_system,
                    start,
                    end,
                    original.len()
                );
            }

            let rev_comp = reverse_complement(&original);

            Ok((start, original, rev_comp))
        }

        NaEdit::Identity => {
            // Identity: same as reference. ref=alt=base at position.
            // We'll encode as a "substitution" where ref==alt.
            let ref_genome = reference.ok_or(HgvsError::ReferenceRequired {
                reason: "identity variants need the reference base",
            })?;

            let base = ref_genome.get_sequence(ucsc_name, start - 1, start)?;
            Ok((start, base.clone(), base))
        }
    }
}

/// Compute the reverse complement of a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

// ── UVID → HGVS conversion ────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct FormatOptions {
    /// If true, attempt to detect duplications and inversions by comparing
    /// the variant against the reference genome. This is more expensive
    /// but produces more informative HGVS notation.
    pub detect_dup_inv: bool,
}

/// Convert a UVID back to HGVS genomic notation.
///
/// # Arguments
///
/// * `uvid` - The UVID to convert.
/// * `reference` - Optional reference genome for dup/inv detection.
/// * `options` - Formatting options.
///
/// # Returns
///
/// The HGVS string, and optionally a list of warnings (e.g., length-mode
/// alleles where exact sequence is unavailable).
pub fn uvid_to_hgvs(
    uvid: &Uvid128,
    reference: Option<&mut dyn ReferenceGenome>,
    options: &FormatOptions,
) -> Result<(String, Vec<String>), HgvsError> {
    let variant = uvid.decode().ok_or(HgvsError::DecodeFailed)?;
    let mut warnings: Vec<String> = Vec::new();

    // Warn about length-mode alleles
    if !variant.ref_is_exact {
        warnings.push(format!(
            "REF allele ({} bp) is in length-mode; using N-repeats as approximation",
            variant.ref_len
        ));
    }
    if !variant.alt_is_exact {
        warnings.push(format!(
            "ALT allele ({} bp) is in length-mode; using N-repeats as approximation",
            variant.alt_len
        ));
    }

    // Get RefSeq accession for this chromosome + assembly
    let chroms = assembly::chromosomes(variant.assembly);
    let chr_idx = variant.chr.0 as usize;
    if chr_idx >= chroms.len() {
        return Err(HgvsError::DecodeFailed);
    }
    let accession = chroms[chr_idx].refseq;

    // Determine coordinate system from chromosome
    let coord_sys = if variant.chr == ChrIndex(24) {
        CoordinateSystem::Mitochondrial
    } else {
        CoordinateSystem::Genomic
    };

    // Convert VCF-style (pos, ref, alt) to HGVS edit
    let hgvs = vcf_to_hgvs_edit(
        &variant,
        accession,
        coord_sys,
        reference,
        options,
        &mut warnings,
    )?;

    Ok((hgvs.to_string(), warnings))
}

/// Convert VCF-style variant fields to an HGVS edit.
fn vcf_to_hgvs_edit(
    variant: &Variant,
    accession: &str,
    coord_sys: CoordinateSystem,
    reference: Option<&mut dyn ReferenceGenome>,
    options: &FormatOptions,
    _warnings: &mut Vec<String>,
) -> Result<HgvsVariant, HgvsError> {
    let ref_seq = &variant.ref_seq;
    let alt_seq = &variant.alt_seq;
    let pos = variant.pos;

    // Determine variant type from VCF alleles
    let ref_len = ref_seq.len();
    let alt_len = alt_seq.len();

    // SNV / MNV: same length, both non-empty
    if ref_len == alt_len && ref_len > 0 {
        if ref_len == 1 {
            // Check for identity
            if ref_seq[0] == alt_seq[0] {
                return Ok(HgvsVariant {
                    accession: accession.to_string(),
                    coordinate_system: coord_sys,
                    start: pos as u64,
                    end: pos as u64,
                    edit: NaEdit::Identity,
                });
            }
            // SNV
            return Ok(HgvsVariant {
                accession: accession.to_string(),
                coordinate_system: coord_sys,
                start: pos as u64,
                end: pos as u64,
                edit: NaEdit::Substitution {
                    reference: ref_seq[0],
                    alternative: alt_seq[0],
                },
            });
        }

        // Optionally detect inversion: alt == reverse complement of ref
        if options.detect_dup_inv && ref_len > 1 {
            let rc = reverse_complement(ref_seq);
            if *alt_seq == rc {
                if ref_len > 12 {
                    _warnings.push(
                        "Inversion allele >12 bp; if stored in length-mode, \
                         sequence is approximate"
                            .to_string(),
                    );
                }
                return Ok(HgvsVariant {
                    accession: accession.to_string(),
                    coordinate_system: coord_sys,
                    start: pos as u64,
                    end: (pos as u64) + (ref_len as u64) - 1,
                    edit: NaEdit::Inversion,
                });
            }
        }

        // MNV: treat as delins (HGVS doesn't have a native MNV notation)
        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start: pos as u64,
            end: (pos as u64) + (ref_len as u64) - 1,
            edit: NaEdit::Delins {
                deleted: Some(ref_seq.clone()),
                inserted: alt_seq.clone(),
            },
        });
    }

    // Pure insertion: ref is anchor-only (1 base), alt starts with anchor
    if ref_len == 1 && alt_len > 1 && alt_seq[0] == ref_seq[0] {
        let inserted = alt_seq[1..].to_vec();

        // Optionally detect duplication
        if options.detect_dup_inv {
            if let Some(ref_genome) = reference {
                let ucsc_name = variant.chr.to_ucsc_name().ok_or(HgvsError::DecodeFailed)?;
                let ins_len = inserted.len() as u32;
                // Check if inserted seq matches the reference at the next position
                // (forward duplication)
                if let Ok(downstream) = ref_genome.get_sequence(ucsc_name, pos, pos + ins_len) {
                    if downstream == inserted {
                        return Ok(HgvsVariant {
                            accession: accession.to_string(),
                            coordinate_system: coord_sys,
                            start: (pos + 1) as u64,
                            end: (pos as u64) + (ins_len as u64),
                            edit: NaEdit::Duplication {
                                sequence: Some(inserted),
                            },
                        });
                    }
                }
            }
        }

        // Plain insertion: between pos and pos+1
        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start: pos as u64,
            end: (pos + 1) as u64,
            edit: NaEdit::Insertion { sequence: inserted },
        });
    }

    // Pure deletion: alt is anchor-only (1 base), ref starts with anchor
    if alt_len == 1 && ref_len > 1 && ref_seq[0] == alt_seq[0] {
        let deleted = ref_seq[1..].to_vec();
        let del_start = pos + 1; // first deleted base (after anchor)
        let del_end = pos + (ref_len as u32) - 1; // last deleted base

        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start: del_start as u64,
            end: del_end as u64,
            edit: NaEdit::Deletion {
                sequence: if variant.ref_is_exact {
                    Some(deleted)
                } else {
                    None
                },
            },
        });
    }

    // Empty alt (deletion without anchor — happens with "." alt in VCF)
    if alt_len == 0 && ref_len > 0 {
        return Ok(HgvsVariant {
            accession: accession.to_string(),
            coordinate_system: coord_sys,
            start: pos as u64,
            end: (pos as u64) + (ref_len as u64) - 1,
            edit: NaEdit::Deletion {
                sequence: if variant.ref_is_exact {
                    Some(ref_seq.clone())
                } else {
                    None
                },
            },
        });
    }

    // Complex: different lengths, not a clean insertion or deletion.
    // Fallback: delins
    Ok(HgvsVariant {
        accession: accession.to_string(),
        coordinate_system: coord_sys,
        start: pos as u64,
        end: if ref_len > 0 {
            (pos as u64) + (ref_len as u64) - 1
        } else {
            pos as u64
        },
        edit: NaEdit::Delins {
            deleted: if variant.ref_is_exact && ref_len > 0 {
                Some(ref_seq.clone())
            } else {
                None
            },
            inserted: alt_seq.clone(),
        },
    })
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::InlineReference;
    use std::collections::HashMap;

    // ── Parser tests ───────────────────────────────────────────────────

    #[test]
    fn test_parse_substitution() {
        let v = parse("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(v.accession, "NC_000001.11");
        assert_eq!(v.coordinate_system, CoordinateSystem::Genomic);
        assert_eq!(v.start, 12345);
        assert_eq!(v.end, 12345);
        assert_eq!(
            v.edit,
            NaEdit::Substitution {
                reference: b'A',
                alternative: b'G',
            }
        );
    }

    #[test]
    fn test_parse_substitution_lowercase() {
        let v = parse("NC_000001.11:g.100a>t").unwrap();
        assert_eq!(
            v.edit,
            NaEdit::Substitution {
                reference: b'A',
                alternative: b'T',
            }
        );
    }

    #[test]
    fn test_parse_deletion_no_seq() {
        let v = parse("NC_000001.11:g.12345_12347del").unwrap();
        assert_eq!(v.start, 12345);
        assert_eq!(v.end, 12347);
        assert_eq!(v.edit, NaEdit::Deletion { sequence: None });
    }

    #[test]
    fn test_parse_deletion_with_seq() {
        let v = parse("NC_000001.11:g.12345_12347delACG").unwrap();
        assert_eq!(v.start, 12345);
        assert_eq!(v.end, 12347);
        assert_eq!(
            v.edit,
            NaEdit::Deletion {
                sequence: Some(b"ACG".to_vec()),
            }
        );
    }

    #[test]
    fn test_parse_single_base_deletion() {
        let v = parse("NC_000001.11:g.12345del").unwrap();
        assert_eq!(v.start, 12345);
        assert_eq!(v.end, 12345);
        assert_eq!(v.edit, NaEdit::Deletion { sequence: None });
    }

    #[test]
    fn test_parse_insertion() {
        let v = parse("NC_000001.11:g.12345_12346insACGT").unwrap();
        assert_eq!(v.start, 12345);
        assert_eq!(v.end, 12346);
        assert_eq!(
            v.edit,
            NaEdit::Insertion {
                sequence: b"ACGT".to_vec(),
            }
        );
    }

    #[test]
    fn test_parse_delins_no_del_seq() {
        let v = parse("NC_000001.11:g.12345_12347delinsACGT").unwrap();
        assert_eq!(v.start, 12345);
        assert_eq!(v.end, 12347);
        assert_eq!(
            v.edit,
            NaEdit::Delins {
                deleted: None,
                inserted: b"ACGT".to_vec(),
            }
        );
    }

    #[test]
    fn test_parse_delins_with_del_seq() {
        let v = parse("NC_000001.11:g.12345_12347delACGinsT").unwrap();
        assert_eq!(v.start, 12345);
        assert_eq!(v.end, 12347);
        assert_eq!(
            v.edit,
            NaEdit::Delins {
                deleted: Some(b"ACG".to_vec()),
                inserted: b"T".to_vec(),
            }
        );
    }

    #[test]
    fn test_parse_duplication_no_seq() {
        let v = parse("NC_000001.11:g.12345_12347dup").unwrap();
        assert_eq!(v.start, 12345);
        assert_eq!(v.end, 12347);
        assert_eq!(v.edit, NaEdit::Duplication { sequence: None });
    }

    #[test]
    fn test_parse_duplication_with_seq() {
        let v = parse("NC_000001.11:g.12345_12347dupACG").unwrap();
        assert_eq!(
            v.edit,
            NaEdit::Duplication {
                sequence: Some(b"ACG".to_vec()),
            }
        );
    }

    #[test]
    fn test_parse_inversion() {
        let v = parse("NC_000001.11:g.12345_12347inv").unwrap();
        assert_eq!(v.start, 12345);
        assert_eq!(v.end, 12347);
        assert_eq!(v.edit, NaEdit::Inversion);
    }

    #[test]
    fn test_parse_identity() {
        let v = parse("NC_000001.11:g.12345=").unwrap();
        assert_eq!(v.start, 12345);
        assert_eq!(v.edit, NaEdit::Identity);
    }

    #[test]
    fn test_parse_mitochondrial() {
        let v = parse("NC_012920.1:m.8993T>G").unwrap();
        assert_eq!(v.coordinate_system, CoordinateSystem::Mitochondrial);
        assert_eq!(v.accession, "NC_012920.1");
        assert_eq!(v.start, 8993);
    }

    #[test]
    fn test_parse_unsupported_coordinate_systems() {
        let result = parse("NM_000001.3:c.100A>G");
        assert!(
            matches!(result, Err(HgvsError::Parse(..))),
            "Expected parse error for NM_ accession, got: {:?}",
            result
        );

        // With NC_ prefix but c. coordinate system
        let result = parse("NC_000001.11:c.100A>G");
        assert!(matches!(
            result,
            Err(HgvsError::UnsupportedCoordinateSystem { system: 'c', .. })
        ));

        let result = parse("NC_000001.11:p.Ala100Gly");
        assert!(matches!(
            result,
            Err(HgvsError::UnsupportedCoordinateSystem { system: 'p', .. })
        ));
    }

    #[test]
    fn test_parse_missing_colon() {
        assert!(parse("NC_000001.11g.12345A>G").is_err());
    }

    #[test]
    fn test_parse_empty_insertion() {
        assert!(parse("NC_000001.11:g.12345_12346ins").is_err());
    }

    // ── Formatting tests ───────────────────────────────────────────────

    #[test]
    fn test_format_substitution() {
        let v = parse("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(v.to_string(), "NC_000001.11:g.12345A>G");
    }

    #[test]
    fn test_format_deletion() {
        let v = parse("NC_000001.11:g.12345_12347del").unwrap();
        assert_eq!(v.to_string(), "NC_000001.11:g.12345_12347del");
    }

    #[test]
    fn test_format_insertion() {
        let v = parse("NC_000001.11:g.12345_12346insACGT").unwrap();
        assert_eq!(v.to_string(), "NC_000001.11:g.12345_12346insACGT");
    }

    #[test]
    fn test_format_delins() {
        let v = parse("NC_000001.11:g.12345_12347delinsACGT").unwrap();
        assert_eq!(v.to_string(), "NC_000001.11:g.12345_12347delinsACGT");
    }

    #[test]
    fn test_format_inversion() {
        let v = parse("NC_000001.11:g.12345_12347inv").unwrap();
        assert_eq!(v.to_string(), "NC_000001.11:g.12345_12347inv");
    }

    // ── Assembly detection tests ───────────────────────────────────────

    #[test]
    fn test_detect_assembly_grch38() {
        let (asm, chr) = detect_assembly("NC_000001.11").unwrap();
        assert_eq!(asm, Assembly::GRCh38);
        assert_eq!(chr, ChrIndex(0));
    }

    #[test]
    fn test_detect_assembly_grch37() {
        let (asm, chr) = detect_assembly("NC_000001.10").unwrap();
        assert_eq!(asm, Assembly::GRCh37);
        assert_eq!(chr, ChrIndex(0));
    }

    #[test]
    fn test_detect_assembly_chrm() {
        // chrM accession is the same in both assemblies; should prefer GRCh38
        let (asm, chr) = detect_assembly("NC_012920.1").unwrap();
        assert_eq!(asm, Assembly::GRCh38);
        assert_eq!(chr, ChrIndex(24));
    }

    #[test]
    fn test_detect_assembly_unknown() {
        assert!(detect_assembly("NC_999999.1").is_err());
    }

    // ── HGVS → UVID conversion tests ──────────────────────────────────

    #[test]
    fn test_hgvs_to_uvid_substitution() {
        // Substitutions are reference-free
        let uvid = hgvs_to_uvid("NC_000001.11:g.12345A>G", None, None).unwrap();
        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.chr, ChrIndex(0));
        assert_eq!(decoded.pos, 12345);
        assert_eq!(decoded.ref_seq, b"A");
        assert_eq!(decoded.alt_seq, b"G");
        assert_eq!(decoded.assembly, Assembly::GRCh38);
    }

    #[test]
    fn test_hgvs_to_uvid_substitution_grch37() {
        let uvid = hgvs_to_uvid("NC_000001.10:g.12345A>G", None, None).unwrap();
        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.assembly, Assembly::GRCh37);
    }

    #[test]
    fn test_hgvs_to_uvid_assembly_mismatch() {
        let result = hgvs_to_uvid("NC_000001.11:g.12345A>G", None, Some(Assembly::GRCh37));
        assert!(matches!(result, Err(HgvsError::AssemblyMismatch { .. })));
    }

    #[test]
    fn test_hgvs_to_uvid_deletion_needs_reference() {
        let result = hgvs_to_uvid("NC_000001.11:g.12345_12347del", None, None);
        assert!(matches!(result, Err(HgvsError::ReferenceRequired { .. })));
    }

    #[test]
    fn test_hgvs_to_uvid_deletion_with_reference() {
        // Build a fake reference for chr1 with enough sequence at position 12344-12347
        // Reference: positions 0-based [12343..12347] = "XACG" where X is the anchor
        let mut seq = vec![b'N'; 12348];
        seq[12343] = b'T'; // anchor at pos 12344
        seq[12344] = b'A'; // pos 12345
        seq[12345] = b'C'; // pos 12346
        seq[12346] = b'G'; // pos 12347

        let mut reference = InlineReference::new(HashMap::new()).with_contig("chr1", &seq);

        let uvid =
            hgvs_to_uvid("NC_000001.11:g.12345_12347del", Some(&mut reference), None).unwrap();

        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.pos, 12344); // anchor position
        assert_eq!(decoded.ref_seq, b"TACG"); // anchor + deleted
        assert_eq!(decoded.alt_seq, b"T"); // anchor only
    }

    #[test]
    fn test_hgvs_to_uvid_insertion_with_reference() {
        let mut seq = vec![b'N'; 12347];
        seq[12344] = b'A'; // anchor at pos 12345

        let mut reference = InlineReference::new(HashMap::new()).with_contig("chr1", &seq);

        let uvid = hgvs_to_uvid(
            "NC_000001.11:g.12345_12346insGGG",
            Some(&mut reference),
            None,
        )
        .unwrap();

        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.pos, 12345);
        assert_eq!(decoded.ref_seq, b"A"); // anchor
        assert_eq!(decoded.alt_seq, b"AGGG"); // anchor + inserted
    }

    #[test]
    fn test_hgvs_to_uvid_delins_with_reference() {
        let mut seq = vec![b'N'; 12348];
        seq[12344] = b'A'; // pos 12345
        seq[12345] = b'C'; // pos 12346
        seq[12346] = b'G'; // pos 12347

        let mut reference = InlineReference::new(HashMap::new()).with_contig("chr1", &seq);

        let uvid = hgvs_to_uvid(
            "NC_000001.11:g.12345_12347delinsT",
            Some(&mut reference),
            None,
        )
        .unwrap();

        let decoded = uvid.decode().unwrap();
        assert_eq!(decoded.pos, 12345);
        assert_eq!(decoded.ref_seq, b"ACG");
        assert_eq!(decoded.alt_seq, b"T");
    }

    // ── UVID → HGVS conversion tests ──────────────────────────────────

    #[test]
    fn test_uvid_to_hgvs_substitution() {
        let uvid = Uvid128::encode(ChrIndex(0), 12345, b"A", b"G", Assembly::GRCh38).unwrap();
        let (hgvs, warnings) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert_eq!(hgvs, "NC_000001.11:g.12345A>G");
        assert!(warnings.is_empty());
    }

    #[test]
    fn test_uvid_to_hgvs_substitution_grch37() {
        let uvid = Uvid128::encode(ChrIndex(0), 12345, b"A", b"G", Assembly::GRCh37).unwrap();
        let (hgvs, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert_eq!(hgvs, "NC_000001.10:g.12345A>G");
    }

    #[test]
    fn test_uvid_to_hgvs_chrx() {
        let uvid = Uvid128::encode(ChrIndex(22), 100, b"C", b"T", Assembly::GRCh38).unwrap();
        let (hgvs, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert_eq!(hgvs, "NC_000023.11:g.100C>T");
    }

    #[test]
    fn test_uvid_to_hgvs_chrm() {
        let uvid = Uvid128::encode(ChrIndex(24), 8993, b"T", b"G", Assembly::GRCh38).unwrap();
        let (hgvs, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        // chrM uses m. coordinate system
        assert_eq!(hgvs, "NC_012920.1:m.8993T>G");
    }

    #[test]
    fn test_uvid_to_hgvs_deletion() {
        // VCF-style deletion: pos=100, ref=TACG, alt=T (anchor T + deleted ACG)
        let uvid = Uvid128::encode(ChrIndex(0), 100, b"TACG", b"T", Assembly::GRCh38).unwrap();
        let (hgvs, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        // Should produce g.101_103delACG (positions after anchor)
        assert_eq!(hgvs, "NC_000001.11:g.101_103delACG");
    }

    #[test]
    fn test_uvid_to_hgvs_insertion() {
        // VCF-style insertion: pos=100, ref=A, alt=AGGG
        let uvid = Uvid128::encode(ChrIndex(0), 100, b"A", b"AGGG", Assembly::GRCh38).unwrap();
        let (hgvs, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        // Should produce g.100_101insGGG
        assert_eq!(hgvs, "NC_000001.11:g.100_101insGGG");
    }

    #[test]
    fn test_uvid_to_hgvs_delins() {
        // VCF-style complex: pos=100, ref=ACG, alt=TT
        let uvid = Uvid128::encode(ChrIndex(0), 100, b"ACG", b"TT", Assembly::GRCh38).unwrap();
        let (hgvs, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert_eq!(hgvs, "NC_000001.11:g.100_102delACGinsTT");
    }

    #[test]
    fn test_uvid_to_hgvs_identity() {
        let uvid = Uvid128::encode(ChrIndex(0), 100, b"A", b"A", Assembly::GRCh38).unwrap();
        let (hgvs, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert_eq!(hgvs, "NC_000001.11:g.100=");
    }

    #[test]
    fn test_uvid_to_hgvs_inversion_detection() {
        // ref=ACG, alt=CGT (reverse complement)
        let uvid = Uvid128::encode(ChrIndex(0), 100, b"ACG", b"CGT", Assembly::GRCh38).unwrap();

        // Without detection: shows as delins
        let (hgvs_no_detect, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert!(
            hgvs_no_detect.contains("del") && hgvs_no_detect.contains("ins"),
            "expected delins notation, got: {}",
            hgvs_no_detect
        );
        assert!(!hgvs_no_detect.contains("inv"));

        // With detection: shows as inv
        let opts = FormatOptions {
            detect_dup_inv: true,
        };
        let (hgvs_detect, _) = uvid_to_hgvs(&uvid, None, &opts).unwrap();
        assert_eq!(hgvs_detect, "NC_000001.11:g.100_102inv");
    }

    #[test]
    fn test_uvid_to_hgvs_length_mode_warning() {
        // 21-base allele triggers length mode
        let ref_seq = b"ACGTACGTACGTACGTACGTA"; // 21 bases
        let uvid = Uvid128::encode(ChrIndex(0), 100, ref_seq, b"T", Assembly::GRCh38).unwrap();
        let (_, warnings) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert!(!warnings.is_empty());
        assert!(warnings[0].contains("length-mode"));
    }

    // ── Round-trip tests ───────────────────────────────────────────────

    #[test]
    fn test_roundtrip_substitution() {
        let input = "NC_000001.11:g.12345A>G";
        let uvid = hgvs_to_uvid(input, None, None).unwrap();
        let (output, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert_eq!(input, output);
    }

    #[test]
    fn test_roundtrip_substitution_chrx() {
        let input = "NC_000023.11:g.100C>T";
        let uvid = hgvs_to_uvid(input, None, None).unwrap();
        let (output, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert_eq!(input, output);
    }

    #[test]
    fn test_roundtrip_mitochondrial() {
        let input = "NC_012920.1:m.8993T>G";
        let uvid = hgvs_to_uvid(input, None, None).unwrap();
        let (output, _) = uvid_to_hgvs(&uvid, None, &FormatOptions::default()).unwrap();
        assert_eq!(input, output);
    }

    // ── Reverse complement tests ───────────────────────────────────────

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ACG"), b"CGT");
        assert_eq!(reverse_complement(b""), b"");
        assert_eq!(reverse_complement(b"A"), b"T");
    }
}
