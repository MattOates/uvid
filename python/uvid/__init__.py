"""UVID - Universal Variant ID for human genetic variation.

Provides compact 128-bit identifiers for genomic variants with
DuckDB-backed collection storage.
"""

from uvid._core import UVID, Collection

__all__ = ["UVID", "Collection"]
__version__ = "0.1.0"
