# Key Concepts

This page explains the core ideas behind UVID encoding.

## Linearized Genome Position

Traditional variant representations use a chromosome name and a position within that chromosome. UVID instead maps all positions into a single **linearized coordinate space** by concatenating chromosomes end-to-end.

Each assembly (GRCh37, GRCh38) has its own offset table. For example, in GRCh38:

- chr1 positions map to offsets 0 -- 248,956,421
- chr2 positions map to offsets 248,956,422 -- 491,149,950
- ...and so on through chrM

This linearized position occupies the top 32 bits of the UVID, which means UVIDs sort in natural genomic order when compared as unsigned 128-bit integers.

!!! note "Bit 127 is not always zero"
    Linearized positions for later chromosomes (e.g. chr14 onwards in GRCh38) exceed 2^31, so the most significant bit of the UVID will be 1. UVIDs must be treated as **unsigned** 128-bit integers for correct sorting.

## Assembly Field

Bits 95-94 encode the genome assembly as a 2-bit value:

| Value | Assembly |
|-------|----------|
| 0 | GRCh37 / hg19 |
| 1 | GRCh38 / hg38 |
| 2-3 | Reserved |

This means UVIDs for GRCh37 and GRCh38 sort into separate groups, with GRCh37 sorting before GRCh38 for the same chromosome and position.

## Allele Encoding

Each allele (REF and ALT) gets a 46-bit field with an independent **mode bit**:

### String Mode (mode = 0)

Used when the allele is 1-20 bases of pure A/C/G/T. Stores the exact sequence.

| Field | Bits | Description |
|-------|------|-------------|
| Mode | 1 | Always 0 |
| Length | 5 | Number of bases (1-20) |
| DNA | 40 | 2-bit encoded bases (A=00, C=01, G=10, T=11) |

String mode alleles sort **before** length mode alleles at the same position, providing a natural ordering where short exact variants come first.

### Length Mode (mode = 1)

Used for alleles longer than 20 bases, or alleles containing non-ACGT characters (e.g. N). Stores the length and a fingerprint.

| Field | Bits | Description |
|-------|------|-------------|
| Mode | 1 | Always 1 |
| Length | 28 | Sequence length in bases (max 268,435,455) |
| Fingerprint | 17 | Rabin fingerprint of the sequence |

Decoding a length-mode allele returns N-repeats as a placeholder. The actual sequence can always be recovered from the reference genome (for REF alleles) or the original VCF.

## Rabin Fingerprint

Length-mode alleles include a 17-bit Rabin fingerprint to reduce collisions between different sequences of the same length at the same locus.

**Polynomial**: x^17 + x^3 + 1 (0x20009 over GF(2))

**Algorithm**: DNA is treated as a 2-bit-per-base bitstream (A=00, C=01, G=10, T=11). For each base, the fingerprint register is shifted left by 2 and the 2-bit base encoding is OR'd in. After each shift, if the register exceeds 17 bits, the polynomial is XOR'd to reduce it.

N bases are treated as A (0b00) -- a zero value requiring no special case.

### Collision resistance

A UVID collision requires two distinct variants to share all of: the same locus (chromosome + position), the same assembly, the same allele mode, and the same allele content encoding. String-mode alleles (up to 20 bases) are stored exactly, so they never collide. Collisions can only occur between length-mode alleles that happen to have the same sequence length **and** the same 17-bit Rabin fingerprint.

The 17-bit fingerprint has 131,072 possible values. For two specific same-length sequences at the same locus, the probability of a fingerprint collision is approximately **1 / 131,072 (~7.6 x 10^-6)**.

Using the birthday-bound approximation, the probability that **any** pair collides among *n* distinct same-length alleles at a single locus is roughly **n^2 / (2 x 131,072)**.

In practice, the number of distinct same-length alleles at any single locus is very small:

| Same-length alleles at one locus | Collision probability |
|----------------------------------|----------------------|
| 2 | ~1 in 131,072 (0.0008%) |
| 5 | ~1 in 10,486 (0.010%) |
| 10 | ~1 in 2,621 (0.038%) |
| 50 | ~1 in 105 (0.95%) |

At genome scale, what matters is whether any locus across an entire dataset produces a collision. The old design (length-only, no fingerprint) had **113 collisions** affecting 265 records across the 4.4 million variant ClinVar dataset. With the 17-bit Rabin fingerprint, the same dataset has **zero collisions** -- the fingerprint provides enough discrimination to separate every real-world case.

!!! success "ClinVar validation"
    4,397,869 ClinVar records processed. Old design: 113 UUID collisions. New design with 17-bit Rabin fingerprint: **zero collisions**.

## UUIDv5

Every UVID can be converted to a standard UUIDv5 for interoperability with systems that expect UUIDs:

```python
uvid = UVID.encode("chr1", 100, "A", "G")
print(uvid.uuid5())  # e.g. a1b2c3d4-e5f6-5789-abcd-ef0123456789
```

The namespace UUID is `2696985c-755c-53de-b6b9-1745af20d0fd`, derived from `uuid5(NAMESPACE_OID, "UVID")`. The name bytes are the raw 128-bit UVID value.

!!! warning "UUIDv5 does not preserve sort order"
    While the UVID itself sorts by genomic position, the UUIDv5 representation does not. Use the raw UVID (hex or integer) when sort order matters.

## Sort Order

UVIDs are designed to sort in a useful order when compared as unsigned 128-bit integers:

1. **Linearized position** (bits 127-96) -- genomic coordinate
2. **Assembly** (bits 95-94) -- GRCh37 before GRCh38
3. **REF mode** (bit 93) -- string mode (exact) before length mode
4. **REF content** (bits 92-47) -- bases or length+fingerprint
5. **ALT mode** (bit 46) -- string mode before length mode
6. **ALT content** (bits 45-0) -- bases or length+fingerprint

This means variants at the same position are grouped together, with short exact alleles sorting before longer ones.
