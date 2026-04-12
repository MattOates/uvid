# UVID

**Compact 128-bit Universal Variant IDs for human genetic variation.**

UVID encodes human genomic variants (SNPs, indels, MNVs) into deterministic 128-bit identifiers that sort by genomic position. Variants are stored in `.uvid` collection files powered by DuckDB, enabling fast region-based queries without a separate database server.

## Features

- **Deterministic encoding** -- the same variant always produces the same 128-bit ID
- **Sort by position** -- UVIDs sort naturally in genomic order as unsigned 128-bit integers
- **Compact** -- 16 bytes per variant, with exact base recovery for alleles up to 20 bp
- **Collision-resistant** -- 17-bit Rabin fingerprint for longer alleles; zero collisions across 4.4 million ClinVar records
- **UUIDv5 compatible** -- every UVID converts to a deterministic UUID for interoperability
- **VCF passthrough** -- stamp UVIDs into VCF ID columns at >200k records/second
- **HGVS support** -- bidirectional conversion between HGVS genomic notation (`g.`/`m.`) and UVIDs
- **DuckDB collections** -- store, search, and query variants by region with no external database

## Quick Example

```python
from uvid import UVID

# Encode a variant
uvid = UVID.encode("chr1", 100, "A", "G", "GRCh38")
print(uvid.to_hex())   # 00000064-40000001-00000000-00000006

# Decode it back
fields = uvid.decode()
# {'chr': '1', 'pos': 100, 'ref': 'A', 'alt': 'G', ...}

# Get a deterministic UUID
print(uvid.uuid5())
```

```bash
# Annotate a VCF
uvid vcf input.vcf output.vcf -a GRCh38

# HGVS to UVID
uvid hgvs-encode "NC_000001.11:g.12345A>G"
```

## Architecture

| Layer | Technology | Role |
|-------|-----------|------|
| Core | Rust | UVID encoding/decoding, VCF parsing (noodles), HGVS notation, DuckDB bulk I/O |
| Bindings | PyO3 + maturin | Zero-copy Python access to Rust core |
| CLI | Typer | Command-line interface wrapping the Rust library |

## Next Steps

- [Installation](getting-started/installation.md) -- install from PyPI or build from source
- [Quick Start](getting-started/quickstart.md) -- encode your first variant
- [Key Concepts](guide/concepts.md) -- understand linearized positions, encoding modes, and Rabin fingerprints
- [Bit Layout](guide/bit-layout.md) -- visual diagrams of the 128-bit structure
