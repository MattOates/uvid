"""Type stubs for the uvid._core native extension module."""

import uuid
from os import PathLike

NAMESPACE_UVID: uuid.UUID
"""The UVID namespace UUID for UUIDv5 generation.

Computed as ``uuid5(NAMESPACE_OID, "UVID")``.
Value: ``2696985c-755c-53de-b6b9-1745af20d0fd``
"""

class AssemblyNotDetectedError(ValueError):
    """Raised when assembly cannot be detected from the VCF header.

    Subclass of ``ValueError`` so it can be caught as either
    ``AssemblyNotDetectedError`` or ``ValueError``.
    """

    ...

class ReferenceNotFoundError(ValueError):
    """Raised when a reference genome file is not found in the data directory.

    Subclass of ``ValueError`` so it can be caught as either
    ``ReferenceNotFoundError`` or ``ValueError``.

    To install reference genomes, run ``uvid setup`` or set
    ``UVID_DATA_DIR`` to a directory containing the reference files.
    """

    ...

def data_dir() -> str | None:
    """Return the platform-specific data directory for UVID reference files.

    Resolution order:

    1. ``UVID_DATA_DIR`` environment variable
    2. Platform default:

       - Linux: ``~/.local/share/uvid/``
       - macOS: ``~/Library/Application Support/uvid/``
       - Windows: ``C:\\Users\\<user>\\AppData\\Roaming\\uvid\\``

    Returns:
        The data directory path, or ``None`` if no platform data
        directory can be determined and ``UVID_DATA_DIR`` is not set.
    """
    ...

def vcf_passthrough(
    input: str | PathLike[str],
    output: str | PathLike[str] | None = None,
    use_uuid: bool = False,
    assembly: str | None = None,
    normalize: bool = False,
) -> int:
    """Process a VCF file, replacing the ID column with UVID identifiers.

    When ``normalize`` is ``True``, variants are normalised using the
    `Tan et al. 2015 <https://doi.org/10.1093/bioinformatics/btv112>`_
    algorithm (left-alignment + parsimonious trimming)
    before UVID encoding, and the output POS/REF/ALT columns reflect
    the normalised representation.  A reference genome file for the
    resolved assembly must be present in the data directory
    (``UVID_DATA_DIR`` or the platform default).

    Args:
        input: Path to input VCF file (.vcf or .vcf.gz).
        output: Path to output file. ``None`` writes plain VCF to stdout.
            If the path ends in ``.vcf.gz``, output is bgzf-compressed.
        use_uuid: When ``True``, emit UUIDv5 representation instead of UVID hex.
        assembly: Assembly override (e.g. ``"GRCh37"``, ``"GRCh38"``).
            ``None`` to auto-detect from the VCF header.
        normalize: When ``True``, normalise variants before encoding.
            Requires a reference genome file in the data directory.

    Returns:
        The number of data records processed.

    Raises:
        AssemblyNotDetectedError: If assembly cannot be detected and no override given.
        OSError: On I/O errors.
        ValueError: On normalisation errors (e.g. reference genome not found).
    """
    ...

def hgvs_to_uvid(
    hgvs: str,
    reference: str | None = None,
    assembly: str | None = None,
) -> UVID:
    """Convert HGVS genomic notation to a UVID.

    Parses an HGVS ``g.`` or ``m.`` expression and encodes it as a
    UVID.  The genome assembly is auto-detected from the RefSeq
    accession version (e.g. ``NC_000001.11`` implies GRCh38).

    A reference genome is required for deletions, insertions, delins,
    duplications, and inversions — only pure substitutions are
    reference-free.

    Args:
        hgvs: HGVS string, e.g. ``"NC_000001.11:g.12345A>G"``.
        reference: Optional path to a ``.2bit`` or ``.fa`` reference
            genome file.  Required for indel-class variants.
        assembly: Optional assembly override for validation (e.g.
            ``"GRCh38"``).  If provided and the accession implies a
            different assembly, a ``ValueError`` is raised.

    Returns:
        A ``UVID`` object.

    Raises:
        ValueError: On parse errors, unknown accessions, assembly
            mismatches, or when a reference genome is required but
            not provided.
    """
    ...

def uvid_to_hgvs(
    uvid: str,
    detect_dup_inv: bool = False,
    reference: str | None = None,
) -> tuple[str, list[str]]:
    """Convert a UVID back to HGVS genomic notation.

    Decodes the UVID's fields and formats them as an HGVS ``g.`` or
    ``m.`` expression.

    Args:
        uvid: Hex string of the UVID (with or without dashes).
        detect_dup_inv: If ``True``, attempt to detect duplications
            and inversions by comparing alleles.  More expensive but
            produces richer notation.  Defaults to ``False``.
        reference: Optional path to a reference genome file.  Required
            when ``detect_dup_inv=True`` for duplication detection.

    Returns:
        A tuple of ``(hgvs_string, warnings)`` where *warnings* is a
        list of strings describing any approximations (e.g.
        length-mode alleles whose exact sequence is unavailable).

    Raises:
        ValueError: If the UVID cannot be decoded or the reference
            genome cannot be opened.
    """
    ...

class UVID:
    """A 128-bit Universal Variant ID encoding a human genomic variant."""

    @staticmethod
    def encode(
        chr: str,
        pos: int,
        ref_seq: str,
        alt_seq: str,
        assembly: str = "GRCh38",
    ) -> UVID:
        """Encode a variant as a UVID.

        Args:
            chr: Chromosome name (e.g. "chr1", "1", "chrX", "X", "chrM").
            pos: 1-based genomic position.
            ref_seq: Reference allele sequence (e.g. "A", "ACGT").
            alt_seq: Alternate allele sequence (e.g. "G", "T", ".").
            assembly: Genome assembly ("GRCh37", "GRCh38", "hg19", "hg38").

        Returns:
            A UVID instance.

        Raises:
            ValueError: If any parameter is invalid.
        """
        ...

    def decode(self) -> dict[str, object]:
        """Decode a UVID back to its component fields.

        Returns:
            A dict with keys: chr (str), pos (int), ref (str), alt (str),
            ref_len (int), alt_len (int), ref_is_exact (bool),
            alt_is_exact (bool), ref_fingerprint (int | None),
            alt_fingerprint (int | None), assembly (str).

            When ``ref_is_exact`` or ``alt_is_exact`` is ``False``, the
            corresponding sequence is returned as N-repeats (the actual
            bases can be recovered from the reference genome).

            ``ref_fingerprint`` and ``alt_fingerprint`` are 17-bit Rabin
            fingerprints of the original sequence, present only for
            length-mode alleles (``None`` for string-mode alleles).

        Raises:
            ValueError: If the UVID data is malformed.
        """
        ...

    def to_hex(self) -> str:
        """Get the hex string representation (format: XXXXXXXX-XXXXXXXX-XXXXXXXX-XXXXXXXX)."""
        ...

    @staticmethod
    def from_hex(hex_str: str) -> UVID:
        """Create a UVID from a hex string (with or without dashes).

        Raises:
            ValueError: If the hex string cannot be parsed.
        """
        ...

    def as_int(self) -> int:
        """Get the raw 128-bit integer value."""
        ...

    @staticmethod
    def from_int(value: int) -> UVID:
        """Create a UVID from a raw 128-bit integer value."""
        ...

    @staticmethod
    def range(
        chr: str,
        start_pos: int,
        end_pos: int,
        assembly: str = "GRCh38",
    ) -> tuple[UVID, UVID]:
        """Compute UVID range bounds for a genomic region.

        Args:
            chr: Chromosome name.
            start_pos: Start position (1-based, inclusive).
            end_pos: End position (1-based, inclusive).
            assembly: Genome assembly ("GRCh37", "GRCh38").

        Returns:
            A (lower, upper) tuple of UVIDs bounding the region.

        Raises:
            ValueError: If the chromosome or position is invalid.
        """
        ...

    def uuid5(self) -> uuid.UUID:
        """Convert this UVID to a deterministic UUIDv5.

        Uses the UVID namespace (derived from OID namespace + "UVID")
        and the raw 128-bit integer bytes as the name.

        Returns:
            A Python uuid.UUID with version=5.
        """
        ...

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __lt__(self, other: UVID) -> bool: ...
    def __le__(self, other: UVID) -> bool: ...
    def __gt__(self, other: UVID) -> bool: ...
    def __ge__(self, other: UVID) -> bool: ...
    def __hash__(self) -> int: ...

class Collection:
    """A .uvid collection file backed by DuckDB."""

    def __init__(self, path: str) -> None:
        """Open or create a .uvid collection file.

        Args:
            path: Path to the .uvid file.

        Raises:
            OSError: If the file cannot be opened or created.
        """
        ...

    def add_vcf(self, vcf_path: str, assembly: str = "GRCh38") -> None:
        """Add a VCF file to the collection.

        Args:
            vcf_path: Path to a VCF file (.vcf or .vcf.gz).
            assembly: Genome assembly ("GRCh37", "GRCh38", "hg19", "hg38").

        Raises:
            OSError: If the VCF file cannot be read or parsed.
            ValueError: If the assembly is invalid.
        """
        ...

    def list_samples(self) -> list[tuple[str, str, str]]:
        """List all samples in the collection.

        Returns:
            A list of (table_name, source_file, sample_name) tuples.

        Raises:
            OSError: If the query fails.
        """
        ...

    def list_sources(self) -> list[str]:
        """List all source files in the collection.

        Returns:
            A list of source file names.

        Raises:
            OSError: If the query fails.
        """
        ...

    def search_region(
        self,
        table_name: str,
        chr: str,
        start_pos: int,
        end_pos: int,
        assembly: str = "GRCh38",
    ) -> list[dict[str, object]]:
        """Search for variants in a genomic region.

        Args:
            table_name: Sample table name (from list_samples).
            chr: Chromosome name.
            start_pos: Start position (1-based, inclusive).
            end_pos: End position (1-based, inclusive).

        Returns:
            A list of dicts with keys: uvid, allele1, allele2, phased,
            dp, gq, qual, filter, multiallelic.

        Raises:
            ValueError: If the chromosome is invalid.
            OSError: If the search fails.
        """
        ...
