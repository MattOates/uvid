# HGVS Notation

UVID supports bidirectional conversion between [HGVS](https://varnomen.hgvs.org/) (Human Genome Variation Society) sequence variant nomenclature and UVIDs.

## Supported Coordinate Systems

Only **genomic** (`g.`) and **mitochondrial** (`m.`) coordinate systems are supported. These map directly to VCF-style chromosomal coordinates used by UVID internally.

Coding (`c.`), non-coding (`n.`), protein (`p.`), and RNA (`r.`) coordinate systems are **not supported** and return a clear error message pointing to the [ferro-hgvs](https://crates.io/crates/ferro-hgvs) crate for transcript-level notation.

## Supported Edit Types

| HGVS Edit | Example | Reference Required |
|-----------|---------|-------------------|
| Substitution | `NC_000001.11:g.12345A>G` | No |
| Identity | `NC_000001.11:g.12345=` | No |
| Deletion | `NC_000001.11:g.12345_12347del` | Yes |
| Insertion | `NC_000001.11:g.12345_12346insACT` | Yes |
| Delins | `NC_000001.11:g.12345_12347delinsACT` | Yes |
| Duplication | `NC_000001.11:g.12345_12347dup` | Yes |
| Inversion | `NC_000001.11:g.12345_12347inv` | Yes |

Substitutions are **reference-free** -- no reference genome is needed. All other edit types require a reference genome to determine anchor bases and deleted sequences for VCF-style encoding.

## Assembly Auto-Detection

The genome assembly is determined from the RefSeq accession version number:

| Accession version | Assembly |
|-------------------|----------|
| `NC_000001.10`, `NC_000002.11`, ... | GRCh37 |
| `NC_000001.11`, `NC_000002.12`, ... | GRCh38 |
| `NC_012920.1` (chrM) | Both (same accession) |

An optional `assembly` parameter can be provided for validation. If the accession implies a different assembly, an error is raised.

## HGVS to UVID

The parser converts HGVS notation to VCF-style coordinates (chromosome, position, REF, ALT) and then encodes them as a UVID.

=== "CLI"

    ```bash
    # Substitution (reference-free)
    uvid hgvs-encode "NC_000001.11:g.12345A>G"

    # Deletion (needs reference)
    uvid hgvs-encode "NC_000001.11:g.12345_12347del" -r hg38.2bit

    # With assembly validation
    uvid hgvs-encode "NC_000001.10:g.12345A>G" -a GRCh37
    ```

=== "Python"

    ```python
    from uvid import hgvs_to_uvid

    # Substitution
    uvid = hgvs_to_uvid("NC_000001.11:g.12345A>G")

    # With reference genome for indels
    uvid = hgvs_to_uvid(
        "NC_000001.11:g.12345_12347del",
        reference="hg38.2bit",
    )

    # With assembly validation
    uvid = hgvs_to_uvid(
        "NC_000001.10:g.12345A>G",
        assembly="GRCh37",
    )
    ```

## UVID to HGVS

Decoding a UVID produces HGVS notation with the RefSeq accession appropriate for the stored assembly.

By default, same-length REF/ALT variants are output as `delins`. Duplication and inversion detection is **opt-in** via the `detect_dup_inv` option, since it requires comparing allele sequences.

=== "CLI"

    ```bash
    # Basic decode
    uvid hgvs-decode 00003039-40000001-00000000-00000006

    # With duplication/inversion detection
    uvid hgvs-decode <uvid-hex> --detect-dup-inv -r hg38.2bit
    ```

=== "Python"

    ```python
    from uvid import uvid_to_hgvs

    hgvs_str, warnings = uvid_to_hgvs("00003039-40000001-00000000-00000006")
    # "NC_000001.11:g.12345A>G"

    # Opt-in dup/inv detection
    hgvs_str, warnings = uvid_to_hgvs(
        uvid_hex,
        detect_dup_inv=True,
        reference="hg38.2bit",
    )
    ```

## Warnings

The conversion functions return a list of warnings in the following cases:

- **Length-mode alleles**: Alleles longer than 20 bp are stored in length-mode in the UVID. The HGVS output uses N-repeats as an approximation since the exact sequence is unavailable from the UVID alone.
- **Inversion length**: When an inversion allele exceeds 12 characters and may have been stored in length-mode, a warning notes that the sequence is approximate.

## Coordinate Mapping

HGVS and VCF use different conventions for representing variants:

| Variant Type | HGVS | VCF |
|-------------|------|-----|
| Substitution | `g.12345A>G` -- position of the changed base | pos=12345, REF=A, ALT=G |
| Insertion | `g.12345_12346insACT` -- flanking positions | pos=12345, REF=\<anchor\>, ALT=\<anchor\>ACT |
| Deletion | `g.12345_12347del` -- range of deleted bases | pos=12344, REF=\<anchor\>XXX, ALT=\<anchor\> |
| Delins | `g.12345_12347delinsACT` -- range + inserted seq | pos=12345, REF=XXX, ALT=ACT |

This is why a reference genome is needed for indels: the anchor base and deleted sequences must be looked up from the reference.
