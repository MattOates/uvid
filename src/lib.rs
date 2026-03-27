//! UVID - Universal Variant ID library with Python bindings.
//!
//! Provides 128-bit compact identifiers for human genetic variation,
//! backed by DuckDB storage and noodles-vcf parsing.

pub mod allele_pack;
pub mod assembly;
pub mod normalize;
pub mod store;
pub mod uvid128;
pub mod vcf;
pub mod vcf_passthrough;

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

use assembly::{classify_chrom, Assembly, ChrIndex};
use uvid128::Uvid128;

// Custom Python exceptions
pyo3::create_exception!(uvid._core, AssemblyNotDetectedError, PyValueError);
pyo3::create_exception!(uvid._core, ReferenceNotFoundError, PyValueError);

/// Python-exposed UVID class.
#[pyclass(name = "UVID", skip_from_py_object)]
#[derive(Clone)]
struct PyUvid {
    inner: Uvid128,
}

#[pymethods]
impl PyUvid {
    /// Create a new UVID by encoding a variant.
    #[staticmethod]
    #[pyo3(signature = (chr, pos, ref_seq, alt_seq, assembly = "GRCh38"))]
    fn encode(chr: &str, pos: u32, ref_seq: &str, alt_seq: &str, assembly: &str) -> PyResult<Self> {
        let asm: Assembly = assembly
            .parse()
            .map_err(|e: String| PyValueError::new_err(e))?;

        let chr_idx = ChrIndex::resolve(chr, classify_chrom(chr), asm)
            .ok_or_else(|| PyValueError::new_err(format!("Invalid chromosome: {}", chr)))?;

        let uvid = Uvid128::encode(chr_idx, pos, ref_seq.as_bytes(), alt_seq.as_bytes(), asm)
            .ok_or_else(|| {
                PyValueError::new_err("Failed to encode UVID: invalid input parameters")
            })?;

        Ok(PyUvid { inner: uvid })
    }

    /// Decode a UVID back to its component fields.
    fn decode<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        let variant = self
            .inner
            .decode()
            .ok_or_else(|| PyValueError::new_err("Failed to decode UVID: malformed data"))?;

        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("chr", variant.chr.to_name().unwrap_or("?"))?;
        dict.set_item("pos", variant.pos)?;
        dict.set_item("ref", String::from_utf8_lossy(&variant.ref_seq).to_string())?;
        dict.set_item("alt", String::from_utf8_lossy(&variant.alt_seq).to_string())?;
        dict.set_item("ref_len", variant.ref_len)?;
        dict.set_item("alt_len", variant.alt_len)?;
        dict.set_item("ref_is_exact", variant.ref_is_exact)?;
        dict.set_item("alt_is_exact", variant.alt_is_exact)?;
        dict.set_item("ref_fingerprint", variant.ref_fingerprint.map(|f| f as u64))?;
        dict.set_item("alt_fingerprint", variant.alt_fingerprint.map(|f| f as u64))?;
        dict.set_item("assembly", variant.assembly.to_string())?;
        Ok(dict)
    }

    /// Get the hex string representation.
    fn to_hex(&self) -> String {
        self.inner.to_hex_string()
    }

    /// Create from a hex string.
    #[staticmethod]
    fn from_hex(hex_str: &str) -> PyResult<Self> {
        let uvid = Uvid128::from_hex_string(hex_str)
            .ok_or_else(|| PyValueError::new_err("Invalid hex string"))?;
        Ok(PyUvid { inner: uvid })
    }

    /// Get the raw integer value as a Python int.
    fn as_int(&self) -> u128 {
        self.inner.as_u128()
    }

    /// Create from a raw integer value.
    #[staticmethod]
    fn from_int(value: u128) -> Self {
        PyUvid {
            inner: Uvid128::from_u128(value),
        }
    }

    /// Compute UVID range bounds for a genomic region.
    ///
    /// Returns (lower, upper) as a tuple of UVID.
    #[staticmethod]
    #[pyo3(signature = (chr, start_pos, end_pos, assembly = "GRCh38"))]
    fn range(
        chr: &str,
        start_pos: u32,
        end_pos: u32,
        assembly: &str,
    ) -> PyResult<(PyUvid, PyUvid)> {
        let asm: Assembly = assembly
            .parse()
            .map_err(|e: String| PyValueError::new_err(e))?;

        let chr_idx = ChrIndex::resolve(chr, classify_chrom(chr), asm)
            .ok_or_else(|| PyValueError::new_err(format!("Invalid chromosome: {}", chr)))?;

        let (lower, upper) = Uvid128::range(chr_idx, start_pos, end_pos, asm).ok_or_else(|| {
            PyValueError::new_err("Failed to compute range: invalid chromosome or position")
        })?;

        Ok((PyUvid { inner: lower }, PyUvid { inner: upper }))
    }

    /// Convert this UVID to a deterministic UUIDv5.
    ///
    /// Uses the UVID namespace (derived from the OID namespace + "UVID")
    /// and the raw 128-bit integer bytes as the name.
    fn uuid5<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let uuid = self.inner.to_uuid5();
        let uuid_mod = py.import("uuid")?;
        uuid_mod.call_method1("UUID", (uuid.to_string(),))
    }

    fn __repr__(&self) -> String {
        format!("UVID({})", self.inner.to_hex_string())
    }

    fn __str__(&self) -> String {
        self.inner.to_hex_string()
    }

    fn __eq__(&self, other: &PyUvid) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    fn __lt__(&self, other: &PyUvid) -> bool {
        self.inner < other.inner
    }

    fn __le__(&self, other: &PyUvid) -> bool {
        self.inner <= other.inner
    }

    fn __gt__(&self, other: &PyUvid) -> bool {
        self.inner > other.inner
    }

    fn __ge__(&self, other: &PyUvid) -> bool {
        self.inner >= other.inner
    }
}

/// Python-exposed Collection class.
#[pyclass(name = "Collection", unsendable)]
struct PyCollection {
    inner: store::UvidStore,
}

#[pymethods]
impl PyCollection {
    /// Open or create a .uvid collection file.
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let store = store::UvidStore::open(&PathBuf::from(path))
            .map_err(|e| PyIOError::new_err(format!("Failed to open store: {}", e)))?;
        Ok(PyCollection { inner: store })
    }

    /// Add a VCF file to the collection.
    #[pyo3(signature = (vcf_path, assembly = "GRCh38"))]
    fn add_vcf(&self, vcf_path: &str, assembly: &str) -> PyResult<()> {
        let path = PathBuf::from(vcf_path);
        let asm: Assembly = assembly
            .parse()
            .map_err(|e: String| PyValueError::new_err(e))?;

        let parsed = crate::vcf::parse_vcf(&path, asm, false)
            .map_err(|e| PyIOError::new_err(format!("Failed to parse VCF: {}", e)))?;

        let prefix = store::UvidStore::table_prefix(&path);

        self.inner
            .add_vcf(&prefix, &parsed)
            .map_err(|e| PyIOError::new_err(format!("Failed to add VCF to store: {}", e)))?;

        Ok(())
    }

    /// List all samples in the collection.
    ///
    /// Returns a list of (table_name, source_file, sample_name) tuples.
    fn list_samples(&self) -> PyResult<Vec<(String, String, String)>> {
        self.inner
            .list_samples()
            .map_err(|e| PyIOError::new_err(format!("Failed to list samples: {}", e)))
    }

    /// List all source files in the collection.
    fn list_sources(&self) -> PyResult<Vec<String>> {
        self.inner
            .list_sources()
            .map_err(|e| PyIOError::new_err(format!("Failed to list sources: {}", e)))
    }

    /// Search for variants in a genomic region.
    #[pyo3(signature = (table_name, chr, start_pos, end_pos, assembly = "GRCh38"))]
    fn search_region<'py>(
        &self,
        py: Python<'py>,
        table_name: &str,
        chr: &str,
        start_pos: u32,
        end_pos: u32,
        assembly: &str,
    ) -> PyResult<Bound<'py, pyo3::types::PyList>> {
        let asm: Assembly = assembly
            .parse()
            .map_err(|e: String| PyValueError::new_err(e))?;

        let chr_idx = ChrIndex::resolve(chr, classify_chrom(chr), asm)
            .ok_or_else(|| PyValueError::new_err(format!("Invalid chromosome: {}", chr)))?;

        let results = self
            .inner
            .search_region(table_name, chr_idx, start_pos, end_pos, asm)
            .map_err(|e| PyIOError::new_err(format!("Search failed: {}", e)))?;

        let list = pyo3::types::PyList::empty(py);
        for r in results {
            let dict = pyo3::types::PyDict::new(py);
            dict.set_item("uvid", r.uvid.to_hex_string())?;
            dict.set_item("allele1", r.allele1)?;
            dict.set_item("allele2", r.allele2)?;
            dict.set_item("phased", r.phased)?;
            dict.set_item("dp", r.dp)?;
            dict.set_item("gq", r.gq)?;
            dict.set_item("qual", r.qual)?;
            dict.set_item("filter", r.filter)?;
            dict.set_item("multiallelic", r.multiallelic)?;
            list.append(dict)?;
        }
        Ok(list)
    }
}

/// The uvid._core Python module.
#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyUvid>()?;
    m.add_class::<PyCollection>()?;

    // Custom exceptions
    m.add(
        "AssemblyNotDetectedError",
        m.py().get_type::<AssemblyNotDetectedError>(),
    )?;
    m.add(
        "ReferenceNotFoundError",
        m.py().get_type::<ReferenceNotFoundError>(),
    )?;

    // Expose the UVID namespace UUID as a Python uuid.UUID constant
    let py = m.py();
    let uuid_mod = py.import("uuid")?;
    let ns = uuid_mod.call_method1("UUID", (uvid128::UVID_NAMESPACE.to_string(),))?;
    m.add("NAMESPACE_UVID", ns)?;

    // Module-level functions
    m.add_function(wrap_pyfunction!(py_vcf_passthrough, m)?)?;
    m.add_function(wrap_pyfunction!(py_data_dir, m)?)?;

    Ok(())
}

/// Process a VCF file, replacing the ID column with UVID identifiers.
///
/// Args:
///     input: Path to input VCF file (.vcf or .vcf.gz).
///     output: Path to output file (None for stdout). If ends in .vcf.gz, bgzf-compressed.
///     use_uuid: If True, emit UUIDv5 instead of UVID hex.
///     assembly: Assembly override ("GRCh37", "GRCh38", etc.). None to auto-detect from header.
///     normalize: If True, normalise variants
///         (Tan et al. 2015, https://doi.org/10.1093/bioinformatics/btv112) before encoding.
///         Requires a reference genome file in the data directory.
///
/// Returns:
///     Number of data records processed.
///
/// Raises:
///     AssemblyNotDetectedError: If assembly cannot be detected and no override given.
///     OSError: On I/O errors.
///     ValueError: On normalization errors (e.g. reference genome not found).
#[pyfunction]
#[pyo3(name = "vcf_passthrough", signature = (input, output = None, use_uuid = false, assembly = None, normalize = false))]
fn py_vcf_passthrough(
    input: PathBuf,
    output: Option<PathBuf>,
    use_uuid: bool,
    assembly: Option<&str>,
    normalize: bool,
) -> PyResult<u64> {
    let asm_override = match assembly {
        Some(s) => {
            let asm: Assembly = s.parse().map_err(|e: String| PyValueError::new_err(e))?;
            Some(asm)
        }
        None => None,
    };

    vcf_passthrough::vcf_passthrough(&input, output.as_deref(), use_uuid, asm_override, normalize)
        .map_err(|e| match e {
            vcf_passthrough::VcfPassthroughError::AssemblyNotDetected => {
                AssemblyNotDetectedError::new_err(e.to_string())
            }
            vcf_passthrough::VcfPassthroughError::Io(io_err) => {
                PyIOError::new_err(format!("{}", io_err))
            }
            vcf_passthrough::VcfPassthroughError::Normalize(ref norm_err) => {
                // Raise a specific ReferenceNotFoundError when the reference
                // genome file is missing, so the CLI can offer to download it.
                if let normalize::NormalizeError::Reference(normalize::ReferenceError::NotFound {
                    ..
                }) = norm_err
                {
                    ReferenceNotFoundError::new_err(format!("{}", norm_err))
                } else {
                    PyValueError::new_err(format!("Normalization error: {}", norm_err))
                }
            }
        })
}

/// Return the platform-specific data directory for UVID reference files.
///
/// Resolution order:
///     1. ``UVID_DATA_DIR`` environment variable
///     2. Platform default (Linux: ``~/.local/share/uvid``, macOS:
///        ``~/Library/Application Support/uvid``, Windows: ``AppData\Roaming\uvid``)
///
/// Returns:
///     The data directory path as a string, or ``None`` if no platform
///     data directory can be determined and ``UVID_DATA_DIR`` is not set.
#[pyfunction]
#[pyo3(name = "data_dir")]
fn py_data_dir() -> Option<String> {
    normalize::data_dir().map(|p| p.to_string_lossy().into_owned())
}
