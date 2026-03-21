//! UVID - Universal Variant ID library with Python bindings.
//!
//! Provides 128-bit compact identifiers for human genetic variation,
//! backed by DuckDB storage and noodles-vcf parsing.

pub mod assembly;
pub mod store;
pub mod twobit;
pub mod uvid128;
pub mod vcf;

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

use assembly::{Assembly, ChrIndex};
use uvid128::Uvid128;

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
        let chr_idx = ChrIndex::from_name(chr)
            .ok_or_else(|| PyValueError::new_err(format!("Invalid chromosome: {}", chr)))?;

        let asm: Assembly = assembly
            .parse()
            .map_err(|e: String| PyValueError::new_err(e))?;

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
        dict.set_item("overflow", variant.overflow)?;
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
    fn range(chr: &str, start_pos: u32, end_pos: u32) -> PyResult<(PyUvid, PyUvid)> {
        let chr_idx = ChrIndex::from_name(chr)
            .ok_or_else(|| PyValueError::new_err(format!("Invalid chromosome: {}", chr)))?;

        let (lower, upper) = Uvid128::range(chr_idx, start_pos, end_pos);

        Ok((PyUvid { inner: lower }, PyUvid { inner: upper }))
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

        let parsed = crate::vcf::parse_vcf(&path, asm)
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
    #[pyo3(signature = (table_name, chr, start_pos, end_pos))]
    fn search_region<'py>(
        &self,
        py: Python<'py>,
        table_name: &str,
        chr: &str,
        start_pos: u32,
        end_pos: u32,
    ) -> PyResult<Bound<'py, pyo3::types::PyList>> {
        let chr_idx = ChrIndex::from_name(chr)
            .ok_or_else(|| PyValueError::new_err(format!("Invalid chromosome: {}", chr)))?;

        let results = self
            .inner
            .search_region(table_name, chr_idx, start_pos, end_pos)
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
    Ok(())
}
