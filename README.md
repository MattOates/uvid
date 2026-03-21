# uvid

Compact 128-bit Universal Variant IDs for human genetic variation, backed by DuckDB storage.

## Overview

`uvid` encodes human genomic variants (SNPs, indels, MNVs) into deterministic 128-bit identifiers that sort by genomic position. Variants are stored in `.uvid` collection files powered by DuckDB, enabling fast region-based queries without a separate database server.

### 128-bit UVID layout (MSB to LSB)

| Bits | Width | Field |
|------|-------|-------|
| 127-123 | 5 | Chromosome (chr1-22=0-21, X=22, Y=23, M=24) |
| 122-91 | 32 | Position (1-based) |
| 90-85 | 6 | REF length (0-63) |
| 84-79 | 6 | ALT length (0-63) |
| 78 | 1 | Overflow flag |
| 77-2 | 76 | Sequence payload (2-bit encoded, 38 bases max) |
| 1-0 | 2 | Assembly (0=GRCh37, 1=GRCh38) |

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

## CLI Usage

```bash
# Encode a variant
uvid encode chr1 100 A G

# Decode a UVID
uvid decode 00000003-20208000-00000000-00000009

# Add a VCF to a collection
uvid add sample.vcf collection.uvid

# List samples
uvid samples collection.uvid

# Search by region
uvid search collection.uvid --sample sample__NA12878 --chr 1 --start 10000 --end 20000

# Collection info
uvid info collection.uvid
```

## Python API

```python
from uvid import UVID, Collection

# Encode/decode
uvid = UVID.encode("chr1", 100, "A", "G", "GRCh38")
print(uvid.to_hex())    # 00000003-20208000-00000000-00000009
fields = uvid.decode()   # {'chr': 'chr1', 'pos': 100, 'ref': 'A', 'alt': 'G', ...}

# Work with .uvid collections
store = Collection("my_variants.uvid")
store.add_vcf("sample.vcf", "GRCh38")
results = store.search_region("sample__NA12878", "chr1", 10000, 20000)
```

## Architecture

- **Rust core**: UVID encoding/decoding, VCF parsing (noodles), DuckDB bulk I/O
- **Python bindings**: PyO3 + maturin
- **CLI**: Typer wrapping the Rust library

## License

Apache-2.0
