//! Reference genome access for variant normalization.
//!
//! Provides a [`ReferenceGenome`] trait and concrete implementations for
//! `.2bit` files ([`TwoBitReference`]) and indexed FASTA files
//! ([`FastaReference`]).  An [`InlineReference`] is provided for unit tests.
//!
//! This module is intentionally free of UVID-specific types.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Seek};
use std::path::Path;

use super::error::ReferenceError;

// ---------------------------------------------------------------------------
// Trait
// ---------------------------------------------------------------------------

/// Read-only access to reference genome sequences.
///
/// Implementations provide random access to a specific genomic region by
/// chromosome name and 0-based half-open coordinates `[start, end)`.
///
/// The returned bytes are uppercase ASCII DNA (`A`, `C`, `G`, `T`, `N`).
pub trait ReferenceGenome {
    /// Fetch the reference sequence for `chrom` in the 0-based half-open
    /// interval `[start, end)`.
    ///
    /// Returns uppercase ASCII DNA bytes.
    fn get_sequence(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<u8>, ReferenceError>;

    /// Return the length of the given chromosome, if known.
    fn chrom_length(&mut self, chrom: &str) -> Result<u32, ReferenceError>;
}

// ---------------------------------------------------------------------------
// TwoBitReference — in-memory .2bit file
// ---------------------------------------------------------------------------

/// Reference genome backed by a `.2bit` file loaded into memory.
///
/// Uses `twobit::TwoBitFile::open_and_read()` for fast random access with
/// no further I/O after the initial load (~800 MB for hg38).
pub struct TwoBitReference {
    inner: twobit::TwoBitMemoryFile,
}

impl TwoBitReference {
    /// Open a `.2bit` file and read it entirely into memory.
    pub fn open(path: &Path) -> Result<Self, ReferenceError> {
        let inner = twobit::TwoBitFile::open_and_read(path)
            .map_err(|e| ReferenceError::Io(std::io::Error::other(e.to_string())))?;
        Ok(TwoBitReference { inner })
    }
}

impl ReferenceGenome for TwoBitReference {
    fn get_sequence(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<u8>, ReferenceError> {
        let seq = self
            .inner
            .read_sequence(chrom, (start as usize)..(end as usize))
            .map_err(|e| {
                let msg = e.to_string();
                if msg.contains("not found") || msg.contains("unknown sequence") {
                    ReferenceError::ChromosomeNotFound(chrom.to_string())
                } else {
                    ReferenceError::Io(std::io::Error::other(msg))
                }
            })?;
        // twobit returns a String; convert to uppercase bytes
        Ok(seq
            .into_bytes()
            .into_iter()
            .map(|b| b.to_ascii_uppercase())
            .collect())
    }

    fn chrom_length(&mut self, chrom: &str) -> Result<u32, ReferenceError> {
        let names = self.inner.chrom_names();
        let sizes = self.inner.chrom_sizes();
        for (name, &size) in names.iter().zip(sizes.iter()) {
            if name == chrom {
                return Ok(size as u32);
            }
        }
        Err(ReferenceError::ChromosomeNotFound(chrom.to_string()))
    }
}

// ---------------------------------------------------------------------------
// FastaReference — indexed FASTA (.fa + .fa.fai)
// ---------------------------------------------------------------------------

/// Reference genome backed by an indexed FASTA file (`.fa` + `.fa.fai`).
///
/// Uses noodles-fasta's `IndexedReader` for seek-based random access.
pub struct FastaReference<R: BufRead + Seek> {
    reader: noodles::fasta::io::IndexedReader<R>,
}

impl FastaReference<BufReader<std::fs::File>> {
    /// Open an indexed FASTA file from disk.
    ///
    /// Expects `path` to point to the `.fa` file; the index is loaded
    /// from `path.with_extension("fa.fai")` (i.e. `ref.fa` → `ref.fa.fai`).
    pub fn open(path: &Path) -> Result<Self, ReferenceError> {
        let fai_path = {
            let mut p = path.as_os_str().to_owned();
            p.push(".fai");
            std::path::PathBuf::from(p)
        };

        let index = noodles::fasta::fai::fs::read(&fai_path).map_err(|e| {
            ReferenceError::Io(std::io::Error::new(
                e.kind(),
                format!("failed to read FASTA index {}: {}", fai_path.display(), e),
            ))
        })?;

        let file = std::fs::File::open(path)?;
        let buf = BufReader::new(file);
        let reader = noodles::fasta::io::IndexedReader::new(buf, index);

        Ok(FastaReference { reader })
    }
}

impl<R: BufRead + Seek> ReferenceGenome for FastaReference<R> {
    fn get_sequence(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<u8>, ReferenceError> {
        use noodles::core::Region;

        // noodles Region uses 1-based inclusive coordinates
        let start_1based =
            noodles::core::Position::try_from((start + 1) as usize).map_err(|_| {
                ReferenceError::OutOfBounds {
                    chrom: chrom.to_string(),
                    requested_start: start,
                    requested_end: end,
                    chrom_length: 0,
                }
            })?;
        let end_1based = noodles::core::Position::try_from(end as usize).map_err(|_| {
            ReferenceError::OutOfBounds {
                chrom: chrom.to_string(),
                requested_start: start,
                requested_end: end,
                chrom_length: 0,
            }
        })?;

        let region = Region::new(chrom, start_1based..=end_1based);

        let record = self.reader.query(&region).map_err(|e| {
            let msg = e.to_string();
            if msg.contains("missing index entry") || msg.contains("invalid reference sequence") {
                ReferenceError::ChromosomeNotFound(chrom.to_string())
            } else {
                ReferenceError::Io(e)
            }
        })?;

        let seq_bytes: Vec<u8> = record
            .sequence()
            .as_ref()
            .iter()
            .map(|&b| b.to_ascii_uppercase())
            .collect();

        Ok(seq_bytes)
    }

    fn chrom_length(&mut self, chrom: &str) -> Result<u32, ReferenceError> {
        let index = self.reader.index();
        for record in index.as_ref() {
            if record.name() == chrom.as_bytes() {
                return Ok(record.length() as u32);
            }
        }
        Err(ReferenceError::ChromosomeNotFound(chrom.to_string()))
    }
}

// ---------------------------------------------------------------------------
// InlineReference — in-memory sequences for unit tests
// ---------------------------------------------------------------------------

/// In-memory reference genome for unit tests.
///
/// Constructed from a `HashMap` of chromosome name → sequence bytes.
/// No file I/O is performed.
pub struct InlineReference {
    sequences: HashMap<String, Vec<u8>>,
}

impl InlineReference {
    /// Create a new inline reference from a map of contig name → sequence.
    pub fn new(sequences: HashMap<String, Vec<u8>>) -> Self {
        InlineReference { sequences }
    }

    /// Convenience builder: add a contig and return `self`.
    pub fn with_contig(mut self, name: &str, seq: &[u8]) -> Self {
        self.sequences.insert(name.to_string(), seq.to_vec());
        self
    }
}

impl ReferenceGenome for InlineReference {
    fn get_sequence(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<u8>, ReferenceError> {
        let seq = self
            .sequences
            .get(chrom)
            .ok_or_else(|| ReferenceError::ChromosomeNotFound(chrom.to_string()))?;

        let start = start as usize;
        let end = end as usize;

        if end > seq.len() {
            return Err(ReferenceError::OutOfBounds {
                chrom: chrom.to_string(),
                requested_start: start as u32,
                requested_end: end as u32,
                chrom_length: seq.len() as u32,
            });
        }

        Ok(seq[start..end].to_vec())
    }

    fn chrom_length(&mut self, chrom: &str) -> Result<u32, ReferenceError> {
        let seq = self
            .sequences
            .get(chrom)
            .ok_or_else(|| ReferenceError::ChromosomeNotFound(chrom.to_string()))?;
        Ok(seq.len() as u32)
    }
}

// ---------------------------------------------------------------------------
// Auto-detect opener
// ---------------------------------------------------------------------------

/// Open a reference genome file, auto-detecting the format by extension.
///
/// - `.2bit` → [`TwoBitReference`] (loaded into memory)
/// - `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz` → [`FastaReference`] (indexed)
///
/// Returns a boxed trait object for uniform handling.
pub fn open_reference(path: &Path) -> Result<Box<dyn ReferenceGenome>, ReferenceError> {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");

    if name.ends_with(".2bit") {
        Ok(Box::new(TwoBitReference::open(path)?))
    } else if name.ends_with(".fa")
        || name.ends_with(".fasta")
        || name.ends_with(".fa.gz")
        || name.ends_with(".fasta.gz")
    {
        Ok(Box::new(FastaReference::open(path)?))
    } else {
        Err(ReferenceError::UnsupportedFormat(path.to_path_buf()))
    }
}
