# uvid

Compact 128-bit Universal Variant IDs for human genetic variation, backed by DuckDB storage.

## Why UVID?

Genomic variant databases typically assign arbitrary integer or string IDs to variants, requiring a round-trip to the database to discover whether a variant already exists before inserting it. UVID eliminates that problem: the ID is computed deterministically from the variant itself, so identical variants always receive identical IDs regardless of where or when they are encoded.

### Benefits

- **Deterministic** -- the same variant always produces the same 128-bit ID, anywhere, without coordinating with a central authority.
- **Streaming-friendly** -- because the ID is known before database interaction, bulk inserts can go straight to upsert in a single transaction without select/update cycles. This dramatically improves throughput for whole-genome batch pipelines.
- **Harmonizing** -- different textual descriptions of the same variant (VCF, HGVS) receive the same ID after normalisation against the reference, eliminating duplicate entries.
- **Sortable** -- UVIDs sort in natural genomic order (chromosome, position, alleles) when compared as unsigned 128-bit integers, so range scans are trivial.
- **Compact** -- 16 bytes per variant. Alleles up to 20 bases are stored exactly; longer alleles keep a 17-bit Rabin fingerprint that eliminated all collisions across 4.4 million ClinVar records.
- **Shard-friendly** -- the ID is known before any database interaction, so it can drive partitioning and sharding strategies in distributed variant stores.
- **UUIDv5 compatible** -- every UVID converts to a deterministic UUIDv5 for interoperability with systems that expect standard UUIDs.

### Limitations

- Only two assemblies are supported: GRCh37/hg19 and GRCh38/hg38 (2 reserved slots remain).
- Alleles longer than 20 bases cannot be stored exactly; decode returns N-repeats as a placeholder (the original sequence is always recoverable from the reference or VCF).
- Focused on human genomics, where the need for good variant infrastructure at clinical scale is most acute.

## Overview

`uvid` encodes human genomic variants (SNPs, indels, MNVs) into deterministic 128-bit identifiers that sort by genomic position. Variants are stored in `.uvid` collection files powered by DuckDB, enabling fast region-based queries without a separate database server.

### 128-bit UVID layout (MSB to LSB)

| Bits | Width | Field |
|------|-------|-------|
| 127-96 | 32 | Linearized genome position |
| 95-94 | 2 | Assembly (0=GRCh37, 1=GRCh38) |
| 93 | 1 | REF mode (0=string, 1=length) |
| 92-47 | 46 | REF payload |
| 46 | 1 | ALT mode (0=string, 1=length) |
| 45-0 | 46 | ALT payload |

**Sort order**: position, assembly, REF mode, REF content, ALT mode, ALT content.

#### Allele encoding (46 bits per allele)

Each allele independently chooses one of two modes:

- **String mode** (mode=0): 5-bit length + 40-bit 2-bit-encoded DNA. Stores up to 20 bases exactly. Sorts before length mode.
- **Length mode** (mode=1): 28-bit length + 17-bit Rabin fingerprint. Used for sequences >20 bases or containing non-ACGT characters. Decode returns N-repeats as a placeholder.

#### Linearized genome position

All chromosomes are concatenated into a single 32-bit coordinate space using assembly-specific offset tables. This replaces the old chromosome+position split and enables natural u128 sorting across the genome.

### UUIDv5

Every UVID can be converted to a deterministic UUIDv5 using the namespace `2696985c-755c-53de-b6b9-1745af20d0fd` (derived from `uuid5(NAMESPACE_OID, "UVID")`).

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
uvid decode 00000064-40000001-00000000-00000006

# Add a VCF to a collection
uvid add collection.uvid sample.vcf

# List samples
uvid samples collection.uvid

# Search by region
uvid search collection.uvid --sample sample__NA12878 --chr chr1 --start 10000 --end 20000

# Annotate a VCF with UVIDs in the ID column
uvid vcf input.vcf output.vcf -a GRCh38

# Collection info
uvid info collection.uvid
```

## Python API

```python
from uvid import UVID, Collection, vcf_passthrough

# Encode/decode
uvid = UVID.encode("chr1", 100, "A", "G", "GRCh38")
print(uvid.to_hex())     # hex string
fields = uvid.decode()    # {'chr': '1', 'pos': 100, 'ref': 'A', 'alt': 'G',
                          #  'ref_len': 1, 'alt_len': 1,
                          #  'ref_is_exact': True, 'alt_is_exact': True,
                          #  'assembly': 'GRCh38'}

# UUIDv5 conversion
print(uvid.uuid5())       # deterministic UUID

# Range queries
lower, upper = UVID.range("chr1", 10000, 20000, "GRCh38")

# Work with .uvid collections
store = Collection("my_variants.uvid")
store.add_vcf("sample.vcf", "GRCh38")
results = store.search_region("sample__NA12878", "chr1", 10000, 20000)

# VCF passthrough (stamps UVIDs into the ID column)
count = vcf_passthrough("input.vcf", "output.vcf", assembly="GRCh38")
```

## Architecture

- **Rust core**: UVID encoding/decoding, VCF parsing (noodles), DuckDB bulk I/O
- **Python bindings**: PyO3 + maturin
- **CLI**: Typer wrapping the Rust library

## License

Apache-2.0
