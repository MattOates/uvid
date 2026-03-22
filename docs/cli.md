# CLI Reference

The `uvid` command-line tool provides access to all core functionality.

## encode

Encode a variant as a UVID.

```bash
uvid encode CHR POS REF ALT [--assembly/-a ASSEMBLY]
```

**Arguments:**

| Argument | Description |
|----------|-------------|
| `CHR` | Chromosome (e.g. `chr1`, `1`, `chrX`, `X`, `chrM`) |
| `POS` | 1-based genomic position |
| `REF` | Reference allele sequence |
| `ALT` | Alternate allele sequence |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--assembly`, `-a` | `GRCh38` | Genome assembly (`GRCh37`, `GRCh38`, `hg19`, `hg38`) |

**Example:**

```bash
$ uvid encode chr1 100 A G
UVID:    00000064-40000001-00000000-00000006
Integer: 527765581483008000
```

## decode

Decode a UVID to its component fields.

```bash
uvid decode UVID_HEX
```

**Arguments:**

| Argument | Description |
|----------|-------------|
| `UVID_HEX` | UVID hex string (with or without dashes) |

**Example:**

```bash
$ uvid decode 00000064-40000001-00000000-00000006
Chromosome: 1
Position:   100
REF:        A (length: 1, exact)
ALT:        G (length: 1, exact)
Assembly:   GRCh38
```

For length-mode alleles, the output includes the Rabin fingerprint:

```bash
$ uvid decode <long-variant-uvid>
Chromosome: 1
Position:   12345
REF:        NNN...N (length: 50, length-only, fingerprint: 0x1a3b5)
ALT:        NNN...N (length: 30, length-only, fingerprint: 0x0f2c8)
Assembly:   GRCh38
```

## vcf

Process a VCF file, replacing the ID column with UVID identifiers.

```bash
uvid vcf INPUT [OUTPUT] [--uuid] [--assembly/-a ASSEMBLY]
```

**Arguments:**

| Argument | Description |
|----------|-------------|
| `INPUT` | Input VCF file (`.vcf` or `.vcf.gz`) |
| `OUTPUT` | Output file; omit for stdout. Use `.vcf.gz` extension for bgzf compression |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--uuid` | `false` | Use UUIDv5 representation instead of UVID hex |
| `--assembly`, `-a` | auto-detect | Override assembly (`GRCh37`, `GRCh38`) |

**Examples:**

```bash
# Basic passthrough
uvid vcf input.vcf output.vcf -a GRCh38

# UUIDv5 format
uvid vcf input.vcf output.vcf --uuid -a GRCh38

# Compressed I/O
uvid vcf input.vcf.gz output.vcf.gz

# Pipe to stdout
uvid vcf input.vcf -a GRCh38 | bgzip > output.vcf.gz
```

## add

Add one or more VCF files to a `.uvid` collection.

```bash
uvid add COLLECTION VCF_PATHS... [--assembly/-a ASSEMBLY]
```

**Arguments:**

| Argument | Description |
|----------|-------------|
| `COLLECTION` | Path to `.uvid` collection file (created if it doesn't exist) |
| `VCF_PATHS` | One or more VCF file paths |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--assembly`, `-a` | `GRCh38` | Genome assembly |

**Example:**

```bash
uvid add my_variants.uvid sample1.vcf sample2.vcf.gz -a GRCh38
```

## samples

List all samples in a `.uvid` collection.

```bash
uvid samples COLLECTION
```

**Example:**

```bash
$ uvid samples my_variants.uvid
Table                                    Source                         Sample
------------------------------------------------------------------------------------------
sample__NA12878                          sample1.vcf                    NA12878
sample__NA12879                          sample2.vcf.gz                 NA12879
```

## search

Search for variants in a genomic region.

```bash
uvid search COLLECTION --sample SAMPLE --chr CHR --start START --end END [--assembly/-a ASSEMBLY]
```

**Options:**

| Option | Description |
|--------|-------------|
| `--sample`, `-s` | Sample table name (from `uvid samples`) |
| `--chr`, `-c` | Chromosome (e.g. `1`, `X`) |
| `--start` | Start position (1-based, inclusive) |
| `--end` | End position (1-based, inclusive) |
| `--assembly`, `-a` | Genome assembly (default: `GRCh38`) |

**Example:**

```bash
$ uvid search my_variants.uvid -s sample__NA12878 -c chr1 --start 10000 --end 20000
Found 3 variant(s):
UVID                                 A1  A2  Ph    DP   GQ Filter
----------------------------------------------------------------------
00002710-40000001-00000000-00000006   0   1   /   150   99 PASS
...
```

## info

Show information about a `.uvid` collection.

```bash
uvid info COLLECTION
```

**Example:**

```bash
$ uvid info my_variants.uvid
Collection: my_variants.uvid
Source files: 2
  - sample1.vcf
  - sample2.vcf.gz
Sample tables: 2
```
