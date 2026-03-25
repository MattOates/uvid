"""Download reference genome files from UCSC."""

from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

# Maps assembly name → (UCSC download URL, local filename).
REFERENCE_URLS: dict[str, tuple[str, str]] = {
    "GRCh38": (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit",
        "GRCh38.2bit",
    ),
    "GRCh37": (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit",
        "GRCh37.2bit",
    ),
}

ASSEMBLIES: list[str] = sorted(REFERENCE_URLS)


def ensure_data_dir(target_dir: Path) -> None:
    """Create the data directory if it doesn't exist."""
    target_dir.mkdir(parents=True, exist_ok=True)


def download_reference(assembly: str, target_dir: Path) -> Path:
    """Download a reference genome file for *assembly* into *target_dir*.

    Uses ``urllib.request`` (stdlib) with a simple line-based progress
    indicator on stderr.  The download is written to a ``.partial`` file
    first, then atomically renamed on success.

    Returns:
        The path to the downloaded file.

    Raises:
        KeyError: If *assembly* is not a recognised name.
        OSError: On network or filesystem errors.
    """
    url, filename = REFERENCE_URLS[assembly]
    dest = target_dir / filename
    partial = dest.with_suffix(dest.suffix + ".partial")

    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as resp:
        total = int(resp.headers.get("Content-Length", 0))
        total_mb = total / (1024 * 1024) if total else 0
        downloaded = 0
        chunk_size = 256 * 1024  # 256 KB

        with partial.open("wb") as fp:
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                fp.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded * 100 // total
                    dl_mb = downloaded / (1024 * 1024)
                    print(
                        f"\rDownloading {filename}... {pct}% ({dl_mb:.0f}/{total_mb:.0f} MB)",
                        end="",
                        file=sys.stderr,
                        flush=True,
                    )
                else:
                    dl_mb = downloaded / (1024 * 1024)
                    print(
                        f"\rDownloading {filename}... {dl_mb:.0f} MB",
                        end="",
                        file=sys.stderr,
                        flush=True,
                    )

        # Newline after the progress line.
        print(file=sys.stderr)

    partial.rename(dest)
    return dest
