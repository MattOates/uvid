"""UVID - Universal Variant ID for human genetic variation.

Provides compact 128-bit identifiers for genomic variants with
DuckDB-backed collection storage.
"""

from uvid._core import (
    NAMESPACE_UVID,
    UVID,
    AssemblyNotDetectedError,
    Collection,
    ReferenceNotFoundError,
    data_dir,
    hgvs_to_uvid,
    uvid_to_hgvs,
    vcf_passthrough,
)

__all__ = [
    "NAMESPACE_UVID",
    "UVID",
    "AssemblyNotDetectedError",
    "Collection",
    "ReferenceNotFoundError",
    "data_dir",
    "hgvs_to_uvid",
    "uvid_to_hgvs",
    "vcf_passthrough",
]
__version__ = "0.1.0"
