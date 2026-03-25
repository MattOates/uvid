"""UVID CLI - command-line interface for managing .uvid collections."""

from pathlib import Path

import typer

from uvid import (
    UVID,
    AssemblyNotDetectedError,
    Collection,
    ReferenceNotFoundError,
    data_dir,
    vcf_passthrough,
)

app = typer.Typer(
    name="uvid",
    help="Universal Variant ID - compact identifiers for human genetic variation.",
    no_args_is_help=True,
)


@app.command()
def setup(
    assembly: list[str] | None = typer.Option(
        None,
        "--assembly",
        "-a",
        help="Assembly to download (GRCh37, GRCh38). Omit for both.",
    ),
) -> None:
    """Download reference genome files for variant normalization.

    Downloads .2bit reference genomes from UCSC into the UVID data
    directory.  Override the data directory with the UVID_DATA_DIR
    environment variable.
    """
    from uvid.download import ASSEMBLIES, REFERENCE_URLS, download_reference, ensure_data_dir

    target = data_dir()
    if target is None:
        typer.echo(
            "Error: Could not determine a data directory.\nSet UVID_DATA_DIR to a writable path.",
            err=True,
        )
        raise typer.Exit(1)

    target_path = Path(target)
    assemblies = list(assembly) if assembly else ASSEMBLIES

    for name in assemblies:
        if name not in REFERENCE_URLS:
            typer.echo(
                f"Error: Unknown assembly '{name}'. Choose from: {', '.join(ASSEMBLIES)}", err=True
            )
            raise typer.Exit(1)

    typer.echo(f"Data directory: {target_path}", err=True)
    ensure_data_dir(target_path)

    for name in assemblies:
        _, filename = REFERENCE_URLS[name]
        dest = target_path / filename
        if dest.exists():
            typer.echo(f"{filename} already exists, skipping.", err=True)
            continue
        typer.echo(f"Downloading {filename} (~800 MB)...", err=True)
        download_reference(name, target_path)
        typer.echo(f"Saved {dest}", err=True)

    typer.echo("Setup complete.", err=True)


@app.command()
def add(
    collection: Path = typer.Argument(..., help="Path to .uvid collection file"),
    vcf_paths: list[Path] = typer.Argument(..., help="Path(s) to VCF file(s) (.vcf or .vcf.gz)"),
    assembly: str = typer.Option(
        "GRCh38", "--assembly", "-a", help="Genome assembly (GRCh37, GRCh38)"
    ),
) -> None:
    """Add one or more VCF files to a .uvid collection."""
    for vcf_path in vcf_paths:
        if not vcf_path.exists():
            typer.echo(f"Error: VCF file not found: {vcf_path}", err=True)
            raise typer.Exit(1)

    store = Collection(str(collection))
    for vcf_path in vcf_paths:
        typer.echo(f"Adding {vcf_path.name} to {collection.name} (assembly: {assembly})...")
        store.add_vcf(str(vcf_path), assembly)
    typer.echo("Done.")


def _handle_reference_not_found(
    exc: ReferenceNotFoundError,
    input: Path,
    output: Path | None,
    uuid: bool,
    assembly: str | None,
    normalize: bool,
) -> int | None:
    """Handle a missing reference genome during ``vcf`` processing.

    If stdin is a TTY, offer to download the reference interactively
    and retry.  Otherwise, print instructions and return ``None``.
    """
    import sys

    from uvid.download import ASSEMBLIES, REFERENCE_URLS, download_reference, ensure_data_dir

    # Try to extract the assembly name from the error message.
    msg = str(exc)
    detected_assembly: str | None = None
    for name in ASSEMBLIES:
        if name in msg:
            detected_assembly = name
            break

    target = data_dir()
    if target is None:
        typer.echo(
            f"Error: {msg}\n\n"
            "Could not determine a data directory.\n"
            "Set UVID_DATA_DIR to a writable path and run:\n"
            "  uvid setup",
            err=True,
        )
        return None

    target_path = Path(target)

    if not sys.stdin.isatty():
        typer.echo(
            f"Error: {msg}\n\n"
            "To install the reference genome, run:\n"
            "  uvid setup\n\n"
            f"Data directory: {target_path}\n"
            "Override with UVID_DATA_DIR environment variable.",
            err=True,
        )
        return None

    # Interactive prompt.
    assembly_label = detected_assembly or "the required assembly"
    typer.echo(
        f"Reference genome for {assembly_label} not found. Download now? (~800 MB) [Y/n] ",
        err=True,
        nl=False,
    )
    answer = sys.stdin.readline().strip().lower()
    if answer and answer not in ("y", "yes"):
        typer.echo(
            "\nTo install later, run:\n"
            "  uvid setup\n\n"
            f"Data directory: {target_path}\n"
            "Override with UVID_DATA_DIR environment variable.",
            err=True,
        )
        return None

    if detected_assembly is None or detected_assembly not in REFERENCE_URLS:
        typer.echo(
            f"Error: Could not determine which assembly to download from: {msg}\n"
            "Run uvid setup --assembly <ASSEMBLY> manually.",
            err=True,
        )
        return None

    # Download and retry.
    ensure_data_dir(target_path)
    download_reference(detected_assembly, target_path)

    try:
        return vcf_passthrough(
            input,
            output,
            use_uuid=uuid,
            assembly=assembly,
            normalize=normalize,
        )
    except ReferenceNotFoundError:
        typer.echo("Error: Reference genome still not found after download.", err=True)
        return None


@app.command()
def vcf(
    input: Path = typer.Argument(..., help="Input VCF file (.vcf or .vcf.gz)"),
    output: Path | None = typer.Argument(
        None,
        help="Output file; omit for stdout. Ends in .vcf.gz for bgzf compression.",
    ),
    uuid: bool = typer.Option(
        False, "--uuid", help="Use UUIDv5 representation instead of UVID hex"
    ),
    assembly: str | None = typer.Option(
        None,
        "--assembly",
        "-a",
        help="Override assembly (GRCh37, GRCh38). Auto-detected from header if omitted.",
    ),
    normalize: bool = typer.Option(
        False,
        "--normalize",
        "-n",
        help=(
            "Normalise variants "
            "(Tan et al. 2015, https://doi.org/10.1093/bioinformatics/btv112) "
            "before encoding. Requires a reference genome in the data directory."
        ),
    ),
) -> None:
    """Process a VCF file, replacing the ID column with UVID identifiers.

    When --normalize is given, variants are left-aligned and trimmed before
    encoding.  The output POS/REF/ALT columns reflect the normalised form.
    A reference genome file for the assembly must be present in the data
    directory (set UVID_DATA_DIR or use the platform default).
    """
    if not input.exists():
        typer.echo(f"Error: VCF file not found: {input}", err=True)
        raise typer.Exit(1)

    try:
        count = vcf_passthrough(
            input,
            output,
            use_uuid=uuid,
            assembly=assembly,
            normalize=normalize,
        )
    except AssemblyNotDetectedError:
        typer.echo(
            "Error: Could not detect genome assembly from VCF header.\n"
            "Use --assembly/-a to specify GRCh37 or GRCh38.",
            err=True,
        )
        raise typer.Exit(1) from None
    except ReferenceNotFoundError as exc:
        result = _handle_reference_not_found(exc, input, output, uuid, assembly, normalize)
        if result is None:
            raise typer.Exit(1) from None
        count = result

    dest = str(output) if output else "stdout"
    typer.echo(f"Processed {count} records → {dest}", err=True)


@app.command()
def samples(
    collection: Path = typer.Argument(..., help="Path to .uvid collection file"),
) -> None:
    """List all samples in a .uvid collection."""
    if not collection.exists():
        typer.echo(f"Error: Collection not found: {collection}", err=True)
        raise typer.Exit(1)

    store = Collection(str(collection))
    sample_list = store.list_samples()

    if not sample_list:
        typer.echo("No samples found.")
        return

    typer.echo(f"{'Table':<40} {'Source':<30} {'Sample':<20}")
    typer.echo("-" * 90)
    for table_name, source_file, sample_name in sample_list:
        typer.echo(f"{table_name:<40} {source_file:<30} {sample_name:<20}")


@app.command()
def info(
    collection: Path = typer.Argument(..., help="Path to .uvid collection file"),
) -> None:
    """Show information about a .uvid collection."""
    if not collection.exists():
        typer.echo(f"Error: Collection not found: {collection}", err=True)
        raise typer.Exit(1)

    store = Collection(str(collection))

    sources = store.list_sources()
    sample_list = store.list_samples()

    typer.echo(f"Collection: {collection.name}")
    typer.echo(f"Source files: {len(sources)}")
    for source in sources:
        typer.echo(f"  - {source}")
    typer.echo(f"Sample tables: {len(sample_list)}")


@app.command()
def search(
    collection: Path = typer.Argument(..., help="Path to .uvid collection file"),
    sample: str = typer.Option(..., "--sample", "-s", help="Sample table name"),
    chr: str = typer.Option(..., "--chr", "-c", help="Chromosome (1-22, X, Y, M)"),
    start: int = typer.Option(..., "--start", help="Start position (1-based, inclusive)"),
    end: int = typer.Option(..., "--end", help="End position (1-based, inclusive)"),
    assembly: str = typer.Option(
        "GRCh38", "--assembly", "-a", help="Genome assembly (GRCh37, GRCh38)"
    ),
) -> None:
    """Search for variants in a genomic region."""
    if not collection.exists():
        typer.echo(f"Error: Collection not found: {collection}", err=True)
        raise typer.Exit(1)

    store = Collection(str(collection))
    results = store.search_region(sample, chr, start, end, assembly)

    if not results:
        typer.echo("No variants found.")
        return

    typer.echo(f"Found {len(results)} variant(s):")
    typer.echo(f"{'UVID':<36} {'A1':>3} {'A2':>3} {'Ph':>3} {'DP':>5} {'GQ':>4} {'Filter':<10}")
    typer.echo("-" * 70)
    for r in results:
        a1 = str(r.get("allele1", ".")) if r.get("allele1") is not None else "."
        a2 = str(r.get("allele2", ".")) if r.get("allele2") is not None else "."
        ph = "|" if r.get("phased") else "/"
        dp = str(r.get("dp", ".")) if r.get("dp") is not None else "."
        gq = str(r.get("gq", ".")) if r.get("gq") is not None else "."
        filt = r.get("filter", ".") or "."
        typer.echo(f"{r['uvid']:<36} {a1:>3} {a2:>3} {ph:>3} {dp:>5} {gq:>4} {filt:<10}")


@app.command()
def encode(
    chr: str = typer.Argument(..., help="Chromosome (1-22, X, Y, M)"),
    pos: int = typer.Argument(..., help="Position (1-based)"),
    ref_seq: str = typer.Argument(..., help="Reference allele sequence"),
    alt_seq: str = typer.Argument(..., help="Alternate allele sequence"),
    assembly: str = typer.Option("GRCh38", "--assembly", "-a", help="Genome assembly"),
) -> None:
    """Encode a variant as a UVID."""
    uvid = UVID.encode(chr, pos, ref_seq, alt_seq, assembly)
    typer.echo(f"UVID:    {uvid.to_hex()}")
    typer.echo(f"Integer: {uvid.as_int()}")


@app.command()
def decode(
    uvid_hex: str = typer.Argument(..., help="UVID hex string (with or without dashes)"),
) -> None:
    """Decode a UVID to its component fields."""
    uvid = UVID.from_hex(uvid_hex)
    fields = uvid.decode()

    typer.echo(f"Chromosome: {fields['chr']}")
    typer.echo(f"Position:   {fields['pos']}")
    ref_mode = "exact" if fields["ref_is_exact"] else "length-only"
    alt_mode = "exact" if fields["alt_is_exact"] else "length-only"
    ref_fp = fields.get("ref_fingerprint")
    alt_fp = fields.get("alt_fingerprint")
    ref_detail = f"length: {fields['ref_len']}, {ref_mode}"
    if ref_fp is not None:
        ref_detail += f", fingerprint: {ref_fp:#07x}"
    alt_detail = f"length: {fields['alt_len']}, {alt_mode}"
    if alt_fp is not None:
        alt_detail += f", fingerprint: {alt_fp:#07x}"
    typer.echo(f"REF:        {fields['ref'] or '.'} ({ref_detail})")
    typer.echo(f"ALT:        {fields['alt'] or '.'} ({alt_detail})")
    typer.echo(f"Assembly:   {fields['assembly']}")


if __name__ == "__main__":
    app()
