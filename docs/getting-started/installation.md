# Installation

UVID requires Python 3.10 or later. The package includes a compiled Rust extension that ships as pre-built wheels for common platforms.

## From PyPI

=== "uv"

    ```bash
    uv pip install uvid
    ```

=== "pip"

    ```bash
    pip install uvid
    ```

## From Source

Building from source requires a Rust toolchain (1.70+). Install Rust via [rustup](https://rustup.rs/) if you don't have it.

```bash
git clone https://github.com/MattOates/uvid.git
cd uvid
uv pip install .
```

This uses [maturin](https://www.maturin.rs/) under the hood to compile the Rust extension and install the Python package.

## As a CLI Tool

If you only need the command-line interface:

=== "uv"

    ```bash
    uv tool install uvid
    ```

=== "pipx"

    ```bash
    pipx install uvid
    ```

This makes the `uvid` command available globally without activating a virtual environment.

## Verify Installation

```bash
uvid --help
```

```python
from uvid import UVID

uvid = UVID.encode("chr1", 100, "A", "G")
print(uvid.to_hex())
```
