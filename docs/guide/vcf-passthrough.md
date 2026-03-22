# VCF Passthrough

The VCF passthrough feature replaces the ID column in a VCF file with UVID identifiers. It is implemented entirely in Rust for performance, processing over 200,000 records per second.

## How It Works

1. Reads the input VCF (plain text or `.vcf.gz` with bgzf compression)
2. Detects the genome assembly from the VCF header (or uses the provided override)
3. For each data record, encodes the variant (CHROM, POS, REF, ALT) as a UVID
4. Replaces the ID column with the UVID hex string (or UUIDv5 if requested)
5. Writes the modified record to the output

Header lines are passed through unchanged. Multi-allelic records generate a UVID from the first ALT allele.

## CLI Usage

```bash
# Basic passthrough
uvid vcf input.vcf output.vcf -a GRCh38

# UUIDv5 format
uvid vcf input.vcf output.vcf -a GRCh38 --uuid

# Auto-detect assembly from header
uvid vcf input.vcf output.vcf

# Compressed I/O
uvid vcf input.vcf.gz output.vcf.gz -a GRCh38

# Stream to stdout
uvid vcf input.vcf -a GRCh38 > output.vcf
uvid vcf input.vcf -a GRCh38 | bgzip > output.vcf.gz
```

## Python API

```python
from uvid import vcf_passthrough

# Basic usage
count = vcf_passthrough("input.vcf", "output.vcf", assembly="GRCh38")
print(f"Processed {count} records")

# UUIDv5 output
count = vcf_passthrough("input.vcf", "output.vcf", use_uuid=True, assembly="GRCh38")

# Auto-detect assembly
count = vcf_passthrough("input.vcf", "output.vcf")

# Stream to stdout (output=None)
count = vcf_passthrough("input.vcf", assembly="GRCh38")

# Compressed output
count = vcf_passthrough("input.vcf", "output.vcf.gz", assembly="GRCh38")
```

## Assembly Detection

When no assembly is specified, the passthrough attempts to detect it from VCF header metadata (e.g. `##reference` or `##contig` lines). If detection fails, an `AssemblyNotDetectedError` is raised:

```python
from uvid import vcf_passthrough, AssemblyNotDetectedError

try:
    vcf_passthrough("input.vcf", "output.vcf")
except AssemblyNotDetectedError:
    # Fall back to explicit assembly
    vcf_passthrough("input.vcf", "output.vcf", assembly="GRCh38")
```

## Performance

Benchmarked on ClinVar (4.4 million records):

| Mode | Throughput |
|------|-----------|
| UVID hex | ~217k records/s |
| UUIDv5 | ~182k records/s |

The bottleneck is UUIDv5 SHA-1 hashing; UVID hex output is faster since it only requires bit manipulation.

## Limitations

- Multi-allelic records use only the first ALT allele for the UVID
- The output VCF retains the original header unchanged (no UVID-specific header lines are added)
- Non-standard contigs (e.g. `NT_113889.1`) are not supported and their records are passed through with the ID column unchanged
