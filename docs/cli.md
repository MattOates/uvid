# CLI Reference

The `uvid` command-line tool provides access to all core functionality.

## setup

Download reference genome files for variant normalization.

```bash
uvid setup [--assembly/-a ASSEMBLY]
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--assembly`, `-a` | both | Assembly to download (`GRCh37`, `GRCh38`). Can be repeated. Omit to download both |

Downloads `.2bit` reference genomes from UCSC into the UVID data directory. The data directory is platform-specific (see [Normalization setup](guide/normalization.md#data-directory)) and can be overridden with the `UVID_DATA_DIR` environment variable.

If a file already exists in the data directory, it is skipped.

**Examples:**

```bash
# Download both GRCh37 and GRCh38
uvid setup

# Download only GRCh38
uvid setup -a GRCh38

# Download to a custom location
UVID_DATA_DIR=/data/references uvid setup -a GRCh38
```

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

## hgvs-encode

Encode an HGVS expression as a UVID.

```bash
uvid hgvs-encode HGVS [--reference/-r PATH] [--assembly/-a ASSEMBLY]
```

Supports genomic (`g.`) and mitochondrial (`m.`) coordinate systems. Substitutions are reference-free; indels (insertions, deletions, delins, duplications, inversions) require a reference genome file.

**Arguments:**

| Argument | Description |
|----------|-------------|
| `HGVS` | HGVS expression (e.g. `NC_000001.11:g.12345A>G`) |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--reference`, `-r` | none | Path to `.2bit` or `.fa` reference genome file. Required for indel-class variants |
| `--assembly`, `-a` | auto-detect | Expected assembly for validation (`GRCh37`, `GRCh38`). Auto-detected from the RefSeq accession version |

The assembly is determined from the RefSeq accession version number: `NC_000001.11` implies GRCh38, `NC_000001.10` implies GRCh37. If `--assembly` is provided, it is validated against the detected assembly and an error is raised on mismatch.

**Examples:**

```bash
# Substitution (no reference needed)
$ uvid hgvs-encode "NC_000001.11:g.12345A>G"
UVID:    00003039-40000001-00000000-00000006
Integer: 207685550972928006

# Deletion (requires reference)
$ uvid hgvs-encode "NC_000001.11:g.12345_12347del" -r hg38.2bit

# Mitochondrial variant
$ uvid hgvs-encode "NC_012920.1:m.8993T>G"

# With assembly validation
$ uvid hgvs-encode "NC_000001.11:g.12345A>G" -a GRCh38
```

## hgvs-decode

Decode a UVID to HGVS genomic notation.

```bash
uvid hgvs-decode UVID_HEX [--detect-dup-inv] [--reference/-r PATH]
```

Produces `g.` or `m.` notation. By default outputs simple `ins`/`del`/`delins` notation; use `--detect-dup-inv` to recognise duplications and inversions.

**Arguments:**

| Argument | Description |
|----------|-------------|
| `UVID_HEX` | UVID hex string (with or without dashes) |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--detect-dup-inv` | `false` | Detect duplications and inversions by comparing alleles. More expensive but produces richer HGVS notation |
| `--reference`, `-r` | none | Path to `.2bit` or `.fa` reference genome file. Needed for duplication detection |

Warnings are printed to stderr when:

- An allele was stored in **length-mode** (>20 bp) and the exact sequence is unavailable -- the output uses N-repeats as an approximation.
- An inversion allele exceeds 12 characters, meaning length-mode encoding may have been used and the sequence is approximate.

**Examples:**

```bash
# Basic decode
$ uvid hgvs-decode 00003039-40000001-00000000-00000006
NC_000001.11:g.12345A>G

# With duplication/inversion detection
$ uvid hgvs-decode <uvid-hex> --detect-dup-inv -r hg38.2bit
NC_000001.11:g.100_102inv
```

## vcf

Process a VCF file, replacing the ID column with UVID identifiers.

```bash
uvid vcf INPUT [OUTPUT] [--uuid] [--normalize/-n] [--assembly/-a ASSEMBLY]
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
| `--normalize`, `-n` | `false` | Normalize variants (left-alignment + trimming) before encoding. Requires a reference genome file; see [Normalization](guide/normalization.md#setup) |
| `--assembly`, `-a` | auto-detect | Override assembly (`GRCh37`, `GRCh38`) |

If `--normalize` is set and the reference genome is missing, `uvid vcf` will
offer to download it interactively when running in a terminal. In non-interactive
contexts (pipes, scripts), it exits with instructions to run `uvid setup`.

**Examples:**

```bash
# Basic passthrough
uvid vcf input.vcf output.vcf -a GRCh38

# With normalization
uvid vcf input.vcf output.vcf --normalize -a GRCh38

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
