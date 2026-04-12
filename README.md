# uvid

[![CI](https://github.com/MattOates/uvid/actions/workflows/ci.yml/badge.svg)](https://github.com/MattOates/uvid/actions/workflows/ci.yml)
[![Docs](https://github.com/MattOates/uvid/actions/workflows/docs.yml/badge.svg)](https://mattoates.github.io/uvid/)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776AB.svg)](https://www.python.org/downloads/)
[![Rust](https://img.shields.io/badge/rust-stable-B7410E.svg)](https://www.rust-lang.org/)

**Compact 128-bit Universal Variant IDs for human genetic variation.**

Encode any human genomic variant -- SNP, indel, MNV -- into a deterministic 128-bit identifier that sorts in natural genomic order. No central authority, no database round-trips, no coordination required.

```python
from uvid import UVID

uvid = UVID.encode("chr1", 100, "A", "G", "GRCh38")
print(uvid)  # 00000064-40000001-00000000-00000006
```

> **[Read the docs](https://mattoates.github.io/uvid/)** for full API reference, normalization guide, and bit-layout details.

---

## Why UVID?

Genomic variant databases typically assign arbitrary integer or string IDs to variants, requiring a round-trip to the database to discover whether a variant already exists before inserting it. UVID eliminates that problem: the ID is computed deterministically from the variant itself, so identical variants always receive identical IDs regardless of where or when they are encoded.

| Property | Detail |
|----------|--------|
| **Deterministic** | Same variant = same ID, anywhere, without coordination |
| **Compact** | 16 bytes per variant (fits in a UUID column) |
| **Sortable** | Natural genomic order when compared as unsigned 128-bit integers |
| **Streaming-friendly** | ID is known before database interaction -- bulk upsert in a single pass |
| **Shard-friendly** | ID-driven partitioning for distributed variant stores |
| **UUIDv5 compatible** | Deterministic SHA-1 mapping for systems that expect standard UUIDs |
| **Sequence-searchable** | Alleles up to 20 bp are stored exactly; longer alleles keep length + a 17-bit Rabin fingerprint |

### 128-bit layout (MSB to LSB)

```
 127      96 95 94 93      47 46      0
 +----------+--+--+---------++---------+
 | position |as|rm| REF     |am| ALT   |
 |  (32)    |(2)(1)| (46)   |(1)| (46) |
 +----------+--+--+---------++---------+
```

| Bits | Width | Field |
|------|-------|-------|
| 127-96 | 32 | Linearized genome position |
| 95-94 | 2 | Assembly (0 = GRCh37, 1 = GRCh38) |
| 93 | 1 | REF mode (0 = string, 1 = length) |
| 92-47 | 46 | REF payload |
| 46 | 1 | ALT mode (0 = string, 1 = length) |
| 45-0 | 46 | ALT payload |

Each allele payload is independently encoded in one of two modes:

- **String mode** (mode=0): 5-bit length + 40-bit 2-bit-encoded DNA. Stores up to 20 bases exactly.
- **Length mode** (mode=1): 28-bit length + 17-bit Rabin fingerprint. Used for sequences >20 bases or containing non-ACGT characters.

### Limitations

- Two assemblies supported: GRCh37/hg19 and GRCh38/hg38 (2 reserved slots remain).
- Alleles longer than 20 bases cannot be decoded exactly; the original sequence is always recoverable from the reference or VCF.
- Focused on human genomics.

---

## Installation

Requires Python 3.10+ and a Rust toolchain (for building from source).

```bash
# From PyPI (once published)
uv pip install uvid

# From source
uv pip install .

# As a CLI tool
uv tool install uvid
```

## Quick Start

### CLI

```bash
# Encode a variant
uvid encode chr1 100 A G

# Decode a UVID
uvid decode 00000064-40000001-00000000-00000006

# Annotate a VCF with UVIDs in the ID column
uvid vcf input.vcf output.vcf -a GRCh38

# Add a VCF to a .uvid collection
uvid add collection.uvid sample.vcf

# Search by region
uvid search collection.uvid --sample sample__NA12878 --chr chr1 --start 10000 --end 20000

# Collection info
uvid info collection.uvid
```

### Python

```python
from uvid import UVID, Collection, vcf_passthrough

# Encode / decode
uvid = UVID.encode("chr1", 100, "A", "G", "GRCh38")
fields = uvid.decode()
# {'chr': '1', 'pos': 100, 'ref': 'A', 'alt': 'G', 'assembly': 'GRCh38', ...}

# UUIDv5 conversion
print(uvid.uuid5())  # deterministic UUID

# Range queries
lower, upper = UVID.range("chr1", 10000, 20000, "GRCh38")

# .uvid collections (DuckDB-backed)
store = Collection("my_variants.uvid")
store.add_vcf("sample.vcf", "GRCh38")
results = store.search_region("sample__NA12878", "chr1", 10000, 20000)

# VCF passthrough -- stamp UVIDs into the ID column
count = vcf_passthrough("input.vcf", "output.vcf", assembly="GRCh38")
```

### Variant Normalization

uvid includes a built-in normalizer based on [Tan et al. 2015](https://doi.org/10.1093/bioinformatics/btv112) (the same algorithm used by bcftools and vt) to ensure consistent IDs across differently-represented variants:

```bash
# Normalize and encode in one pass
uvid vcf input.vcf output.vcf -a GRCh38 --normalize
```

See the [normalization guide](https://mattoates.github.io/uvid/guide/normalization/) for details on reference genome setup.

---

## Architecture

```
                   Python (typer CLI / library API)
                          |
                       PyO3 FFI
                          |
  +-----+--------+-------+--------+-----------+
  |     |        |       |        |           |
  | uvid128  assembly  vcf  normalize  store  |
  | (encode/ (chr     (noodles (Tan et  (DuckDB|
  |  decode)  offsets)  parser) al 2015) I/O)  |
  +--------------------------------------------+
                    Rust core
```

- **Rust core**: UVID encoding/decoding, VCF parsing via [noodles](https://github.com/zaeleus/noodles), variant normalization, DuckDB bulk I/O
- **Python bindings**: [PyO3](https://pyo3.rs) + [maturin](https://www.maturin.rs)
- **CLI**: [Typer](https://typer.tiangolo.com) wrapping the native library

## License

Apache-2.0
