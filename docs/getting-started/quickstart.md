# Quick Start

This guide walks through the core operations: encoding variants, decoding UVIDs, annotating VCFs, and working with collections.

## Encode a Variant

```python
from uvid import UVID

uvid = UVID.encode("chr1", 100, "A", "G", "GRCh38")
print(uvid.to_hex())     # 00000064-40000001-00000000-00000006
print(uvid.as_int())      # raw 128-bit integer
print(uvid.uuid5())       # deterministic UUIDv5
```

All parameters are required except `assembly`, which defaults to `"GRCh38"`. Chromosome names accept both `"chr1"` and `"1"` formats.

## Decode a UVID

```python
uvid = UVID.from_hex("00000064-40000001-00000000-00000006")
fields = uvid.decode()
```

The returned dict contains:

| Key | Type | Description |
|-----|------|-------------|
| `chr` | `str` | Chromosome number (e.g. `"1"`, `"X"`) |
| `pos` | `int` | 1-based genomic position |
| `ref` | `str` | Reference allele sequence |
| `alt` | `str` | Alternate allele sequence |
| `ref_len` | `int` | Reference allele length in bases |
| `alt_len` | `int` | Alternate allele length in bases |
| `ref_is_exact` | `bool` | `True` if REF was stored as exact bases |
| `alt_is_exact` | `bool` | `True` if ALT was stored as exact bases |
| `ref_fingerprint` | `int \| None` | 17-bit Rabin fingerprint (length-mode only) |
| `alt_fingerprint` | `int \| None` | 17-bit Rabin fingerprint (length-mode only) |
| `assembly` | `str` | Genome assembly name |

When `ref_is_exact` or `alt_is_exact` is `False`, the sequence is returned as N-repeats. The actual bases can be recovered from the reference genome.

## Annotate a VCF

The fastest way to add UVIDs to an existing VCF:

=== "CLI"

    ```bash
    # UVID hex in the ID column
    uvid vcf input.vcf output.vcf -a GRCh38

    # UUIDv5 in the ID column
    uvid vcf input.vcf output.vcf -a GRCh38 --uuid

    # Stream to stdout (pipe-friendly)
    uvid vcf input.vcf -a GRCh38 | bgzip > output.vcf.gz
    ```

=== "Python"

    ```python
    from uvid import vcf_passthrough

    count = vcf_passthrough("input.vcf", "output.vcf", assembly="GRCh38")
    print(f"Processed {count} records")

    # UUIDv5 mode
    count = vcf_passthrough("input.vcf", "output.vcf", use_uuid=True)

    # Auto-detect assembly from VCF header
    count = vcf_passthrough("input.vcf", "output.vcf")
    ```

The passthrough processes >200k records/second and supports `.vcf.gz` input and output.

## Work with Collections

Collections are `.uvid` files backed by DuckDB for fast region-based queries.

```python
from uvid import Collection

# Create or open a collection
store = Collection("my_variants.uvid")

# Add a VCF
store.add_vcf("sample.vcf", "GRCh38")

# List samples
for table_name, source_file, sample_name in store.list_samples():
    print(f"{table_name}: {sample_name} from {source_file}")

# Search a region
results = store.search_region(
    "sample__NA12878", "chr1", 10000, 20000, "GRCh38"
)
for r in results:
    print(r["uvid"], r["allele1"], r["allele2"])
```

## Range Queries

Compute UVID bounds for a genomic region (useful for database range scans):

```python
lower, upper = UVID.range("chr1", 10000, 20000, "GRCh38")
# All UVIDs in [lower, upper] fall within chr1:10000-20000
```

## Next Steps

- [Key Concepts](../guide/concepts.md) -- understand linearized positions and encoding modes
- [Bit Layout](../guide/bit-layout.md) -- visual packet diagrams of the 128-bit structure
- [CLI Reference](../cli.md) -- full command-line usage
- [Python API](../reference/python-api.md) -- complete API documentation
