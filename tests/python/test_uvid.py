"""Tests for the UVID Python bindings and CLI."""

from __future__ import annotations

import uuid
from pathlib import Path

import pytest
from typer.testing import CliRunner

from uvid import (
    NAMESPACE_UVID,
    UVID,
    AssemblyNotDetectedError,
    Collection,
    vcf_passthrough,
)
from uvid.cli import app

# ──────────────────────────────────────────────────────
# UVID class tests
# ──────────────────────────────────────────────────────


class TestUVIDEncode:
    """Test UVID encoding."""

    def test_encode_snp(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        assert uvid is not None
        assert str(uvid)  # hex string

    def test_encode_snp_explicit_assembly(self):
        uvid = UVID.encode("chr1", 100, "A", "G", "GRCh38")
        assert uvid is not None

    def test_encode_grch37(self):
        uvid37 = UVID.encode("chr1", 100, "A", "G", "GRCh37")
        uvid38 = UVID.encode("chr1", 100, "A", "G", "GRCh38")
        assert uvid37 != uvid38

    def test_encode_deletion(self):
        uvid = UVID.encode("chr1", 100, "ACGT", "A")
        fields = uvid.decode()
        assert fields["ref"] == "ACGT"
        assert fields["alt"] == "A"

    def test_encode_insertion(self):
        uvid = UVID.encode("chr1", 100, "A", "ACGTACGT")
        fields = uvid.decode()
        assert fields["ref"] == "A"
        assert fields["alt"] == "ACGTACGT"

    def test_encode_mnv(self):
        uvid = UVID.encode("chr1", 100, "ACG", "TGA")
        fields = uvid.decode()
        assert fields["ref"] == "ACG"
        assert fields["alt"] == "TGA"

    def test_encode_all_chromosomes(self):
        for i in range(1, 23):
            uvid = UVID.encode(f"chr{i}", 1, "A", "G")
            assert uvid.decode()["chr"] == str(i)

        uvid_x = UVID.encode("chrX", 1, "A", "G")
        assert uvid_x.decode()["chr"] == "X"

        uvid_y = UVID.encode("chrY", 1, "A", "G")
        assert uvid_y.decode()["chr"] == "Y"

        uvid_m = UVID.encode("chrM", 1, "A", "G")
        assert uvid_m.decode()["chr"] == "M"

    def test_encode_without_chr_prefix(self):
        uvid = UVID.encode("1", 100, "A", "G")
        assert uvid.decode()["chr"] == "1"

    def test_encode_invalid_chromosome(self):
        with pytest.raises(ValueError, match="Invalid chromosome"):
            UVID.encode("chr99", 100, "A", "G")

    def test_encode_n_bases_use_length_mode(self):
        """N-containing sequences no longer fail — they use length mode."""
        uvid = UVID.encode("chr1", 100, "N", "G")
        fields = uvid.decode()
        assert fields["ref_is_exact"] is False
        assert fields["alt_is_exact"] is True
        # Length mode has fingerprint, string mode does not
        assert fields["ref_fingerprint"] is not None
        assert fields["alt_fingerprint"] is None

    def test_encode_position_zero(self):
        with pytest.raises(ValueError):
            UVID.encode("chr1", 0, "A", "G")

    def test_encode_hg19_alias(self):
        """hg19 is a valid alias for GRCh37."""
        uvid = UVID.encode("chr1", 100, "A", "G", "hg19")
        assert uvid.decode()["assembly"] == "GRCh37"

    def test_encode_invalid_assembly(self):
        with pytest.raises(ValueError):
            UVID.encode("chr1", 100, "A", "G", "nonsense_assembly")


class TestUVIDDecode:
    """Test UVID decoding."""

    def test_decode_roundtrip(self):
        uvid = UVID.encode("chr3", 12345, "GCTA", "A")
        fields = uvid.decode()
        assert fields["chr"] == "3"
        assert fields["pos"] == 12345
        assert fields["ref"] == "GCTA"
        assert fields["alt"] == "A"
        assert fields["ref_len"] == 4
        assert fields["alt_len"] == 1
        assert fields["ref_is_exact"] is True
        assert fields["alt_is_exact"] is True
        assert fields["ref_fingerprint"] is None
        assert fields["alt_fingerprint"] is None
        assert fields["assembly"] == "GRCh38"

    def test_decode_length_mode(self):
        """Sequences > 20 bases use length mode (not exact)."""
        uvid = UVID.encode("chr1", 1, "A" * 21, "C" * 21)
        fields = uvid.decode()
        assert fields["ref_is_exact"] is False
        assert fields["alt_is_exact"] is False
        assert fields["ref_len"] == 21
        assert fields["alt_len"] == 21
        # Length mode returns N-repeats as placeholder
        assert fields["ref"] == "N" * 21
        assert fields["alt"] == "N" * 21
        # Length mode has fingerprints
        assert fields["ref_fingerprint"] is not None
        assert fields["alt_fingerprint"] is not None
        assert isinstance(fields["ref_fingerprint"], int)
        assert isinstance(fields["alt_fingerprint"], int)
        # Fingerprints should be < 2^17 = 131072
        assert 0 <= fields["ref_fingerprint"] < 131072
        assert 0 <= fields["alt_fingerprint"] < 131072

    def test_decode_max_string_mode(self):
        """Exactly 20 bases per allele stays in string mode (exact)."""
        uvid = UVID.encode("chr1", 1, "A" * 20, "C" * 20)
        fields = uvid.decode()
        assert fields["ref_is_exact"] is True
        assert fields["alt_is_exact"] is True
        assert fields["ref_fingerprint"] is None
        assert fields["alt_fingerprint"] is None
        assert fields["ref"] == "A" * 20
        assert fields["alt"] == "C" * 20

    def test_decode_fingerprint_discriminates_sequences(self):
        """Same-length sequences at same locus should produce different fingerprints."""
        uvid_a = UVID.encode("chr1", 1, "A", "A" * 21 + "C")  # 22 bases ending in C
        uvid_b = UVID.encode("chr1", 1, "A", "A" * 21 + "G")  # 22 bases ending in G
        fields_a = uvid_a.decode()
        fields_b = uvid_b.decode()
        assert fields_a["alt_len"] == fields_b["alt_len"] == 22
        # Fingerprints should differ
        assert fields_a["alt_fingerprint"] != fields_b["alt_fingerprint"]
        # And therefore the UVIDs themselves should differ
        assert uvid_a != uvid_b


class TestUVIDHex:
    """Test UVID hex string conversion."""

    def test_to_hex(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        hex_str = uvid.to_hex()
        assert isinstance(hex_str, str)
        # Format: 8-8-8-8 hex digits with dashes
        parts = hex_str.split("-")
        assert len(parts) == 4
        assert all(len(p) == 8 for p in parts)

    def test_from_hex_roundtrip(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        hex_str = uvid.to_hex()
        recovered = UVID.from_hex(hex_str)
        assert uvid == recovered

    def test_from_hex_no_dashes(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        hex_str = uvid.to_hex().replace("-", "")
        recovered = UVID.from_hex(hex_str)
        assert uvid == recovered

    def test_from_hex_lenient(self):
        """from_hex strips non-hex chars and parses what remains."""
        # This is lenient by design — non-hex chars are filtered out
        uvid = UVID.from_hex("not-a-valid-hex-string!")
        assert uvid is not None

    def test_from_hex_too_long(self):
        """More than 32 hex digits should fail."""
        with pytest.raises(ValueError, match="Invalid hex string"):
            UVID.from_hex("a" * 33)


class TestUVIDInt:
    """Test UVID integer conversion."""

    def test_as_int(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        value = uvid.as_int()
        assert isinstance(value, int)
        assert value > 0

    def test_from_int_roundtrip(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        value = uvid.as_int()
        recovered = UVID.from_int(value)
        assert uvid == recovered


class TestUVIDSorting:
    """Test UVID sorting properties."""

    def test_chromosome_sorting(self):
        """chr1 < chr2 < ... < chrX < chrY < chrM."""
        uvid1 = UVID.encode("chr1", 100, "A", "G")
        uvid2 = UVID.encode("chr2", 100, "A", "G")
        assert uvid1 < uvid2

    def test_position_sorting(self):
        """Lower position sorts before higher."""
        uvid1 = UVID.encode("chr1", 100, "A", "G")
        uvid2 = UVID.encode("chr1", 200, "A", "G")
        assert uvid1 < uvid2

    def test_assembly_clustering(self):
        """Same-locus variants from different assemblies differ only in assembly bits."""
        uvid37 = UVID.encode("chr1", 100, "A", "G", "GRCh37")
        uvid38 = UVID.encode("chr1", 100, "A", "G", "GRCh38")
        # chr1 has offset 0 in both assemblies, so linearized pos is the same.
        # Assembly is at bits 95-94, so they differ by exactly 1 << 94.
        assert uvid37 != uvid38
        # GRCh37 = assembly code 0, GRCh38 = assembly code 1
        # So GRCh37 UVID < GRCh38 UVID for same locus
        assert uvid37 < uvid38

    def test_list_sort(self):
        """UVIDs should sort in genomic order."""
        uvids = [
            UVID.encode("chr2", 100, "A", "G"),
            UVID.encode("chr1", 200, "A", "G"),
            UVID.encode("chr1", 100, "A", "G"),
        ]
        sorted_uvids = sorted(uvids)
        assert sorted_uvids[0].decode()["chr"] == "1"
        assert sorted_uvids[0].decode()["pos"] == 100
        assert sorted_uvids[1].decode()["chr"] == "1"
        assert sorted_uvids[1].decode()["pos"] == 200
        assert sorted_uvids[2].decode()["chr"] == "2"


class TestUVIDDunderMethods:
    """Test UVID magic methods."""

    def test_repr(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        assert repr(uvid).startswith("UVID(")

    def test_str(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        assert str(uvid) == uvid.to_hex()

    def test_eq(self):
        a = UVID.encode("chr1", 100, "A", "G")
        b = UVID.encode("chr1", 100, "A", "G")
        assert a == b

    def test_ne(self):
        a = UVID.encode("chr1", 100, "A", "G")
        b = UVID.encode("chr1", 100, "A", "T")
        assert a != b

    def test_hash_consistent(self):
        a = UVID.encode("chr1", 100, "A", "G")
        b = UVID.encode("chr1", 100, "A", "G")
        assert hash(a) == hash(b)

    def test_hash_in_set(self):
        a = UVID.encode("chr1", 100, "A", "G")
        b = UVID.encode("chr1", 100, "A", "G")
        c = UVID.encode("chr1", 100, "A", "T")
        s = {a, b, c}
        assert len(s) == 2

    def test_comparison_operators(self):
        a = UVID.encode("chr1", 100, "A", "G")
        b = UVID.encode("chr1", 200, "A", "G")
        assert a < b
        assert a <= b
        assert b > a
        assert b >= a
        assert a <= a
        assert a >= a


class TestUVIDRange:
    """Test UVID range query bounds."""

    def test_range_basic(self):
        lower, upper = UVID.range("chr1", 1000, 2000)
        assert lower < upper

    def test_range_contains(self):
        uvid = UVID.encode("chr3", 1500, "A", "G")
        lower, upper = UVID.range("chr3", 1000, 2000)
        assert lower <= uvid <= upper

    def test_range_excludes(self):
        uvid = UVID.encode("chr3", 500, "A", "G")
        lower, upper = UVID.range("chr3", 1000, 2000)
        assert uvid < lower

    def test_range_invalid_chr(self):
        with pytest.raises(ValueError):
            UVID.range("chr99", 1, 100)


class TestUVIDDeterminism:
    """Test UVID encoding is deterministic."""

    def test_same_input_same_output(self):
        a = UVID.encode("chr1", 100, "A", "G")
        b = UVID.encode("chr1", 100, "A", "G")
        assert a == b
        assert a.as_int() == b.as_int()
        assert a.to_hex() == b.to_hex()

    def test_different_input_different_output(self):
        a = UVID.encode("chr1", 100, "A", "G")
        b = UVID.encode("chr1", 100, "A", "T")
        assert a != b


class TestUVIDUuid5:
    """Test UVID to UUIDv5 conversion."""

    def test_uuid5_returns_uuid(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        result = uvid.uuid5()
        assert isinstance(result, uuid.UUID)

    def test_uuid5_version(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        result = uvid.uuid5()
        assert result.version == 5

    def test_uuid5_rfc4122_variant(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        result = uvid.uuid5()
        # RFC 4122 variant has the top bits of clock_seq_hi_and_reserved == 0b10
        assert result.variant == "specified in RFC 4122"

    def test_uuid5_deterministic(self):
        uvid = UVID.encode("chr1", 100, "A", "G")
        a = uvid.uuid5()
        b = uvid.uuid5()
        assert a == b

    def test_uuid5_different_variants(self):
        a = UVID.encode("chr1", 100, "A", "G")
        b = UVID.encode("chr1", 100, "A", "T")
        assert a.uuid5() != b.uuid5()

    def test_uuid5_matches_python_stdlib(self):
        """UUIDv5 from Rust should match manual SHA-1 construction with same inputs."""
        import hashlib

        uvid = UVID.encode("chr1", 100, "A", "G")
        rust_uuid = uvid.uuid5()
        # Replicate UUIDv5 manually: SHA-1(namespace_bytes + name_bytes)
        # Python's uuid.uuid5() encodes the name as UTF-8 which mangles
        # raw bytes, so we compute it directly per RFC 4122 §4.3.
        raw_bytes = uvid.as_int().to_bytes(16, byteorder="big")
        h = hashlib.sha1(NAMESPACE_UVID.bytes + raw_bytes).digest()[:16]
        b = bytearray(h)
        b[6] = (b[6] & 0x0F) | 0x50  # version 5
        b[8] = (b[8] & 0x3F) | 0x80  # RFC 4122 variant
        python_uuid = uuid.UUID(bytes=bytes(b))
        assert rust_uuid == python_uuid


class TestNamespaceUVID:
    """Test the NAMESPACE_UVID module constant."""

    def test_is_uuid(self):
        assert isinstance(NAMESPACE_UVID, uuid.UUID)

    def test_value(self):
        expected = uuid.uuid5(uuid.NAMESPACE_OID, "UVID")
        assert expected == NAMESPACE_UVID

    def test_version(self):
        assert NAMESPACE_UVID.version == 5

    def test_string_value(self):
        assert str(NAMESPACE_UVID) == "2696985c-755c-53de-b6b9-1745af20d0fd"


# ──────────────────────────────────────────────────────
# Collection class tests
# ──────────────────────────────────────────────────────

TEST_VCF = Path(__file__).parent.parent.parent / "test_data" / "sample.vcf"


class TestCollectionCreate:
    """Test Collection creation and basic operations."""

    def test_create_new(self, tmp_path):
        path = tmp_path / "test.uvid"
        Collection(str(path))  # side effect: creates the file
        assert path.exists()

    def test_open_existing(self, tmp_path):
        path = tmp_path / "test.uvid"
        Collection(str(path))
        # Open again
        store = Collection(str(path))
        samples = store.list_samples()
        assert isinstance(samples, list)

    def test_empty_collection(self, tmp_path):
        path = tmp_path / "test.uvid"
        store = Collection(str(path))
        assert store.list_samples() == []
        assert store.list_sources() == []


class TestCollectionAddVcf:
    """Test adding VCF files to a collection."""

    @pytest.fixture
    def collection(self, tmp_path):
        path = tmp_path / "test.uvid"
        return Collection(str(path))

    def test_add_vcf(self, collection):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        collection.add_vcf(str(TEST_VCF), "GRCh38")

    def test_add_vcf_lists_samples(self, collection):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        collection.add_vcf(str(TEST_VCF), "GRCh38")
        samples = collection.list_samples()
        assert len(samples) > 0
        # Each sample is (table_name, source_file, sample_name)
        assert len(samples[0]) == 3

    def test_add_vcf_lists_sources(self, collection):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        collection.add_vcf(str(TEST_VCF), "GRCh38")
        sources = collection.list_sources()
        assert len(sources) == 1
        assert sources[0] == "sample"

    def test_add_vcf_invalid_path(self, collection):
        with pytest.raises(OSError):
            collection.add_vcf("/nonexistent/file.vcf", "GRCh38")

    def test_add_vcf_invalid_assembly(self, collection):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        with pytest.raises(ValueError):
            collection.add_vcf(str(TEST_VCF), "nonsense_assembly")


class TestCollectionSearch:
    """Test searching a collection."""

    @pytest.fixture
    def loaded_collection(self, tmp_path):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        path = tmp_path / "test.uvid"
        store = Collection(str(path))
        store.add_vcf(str(TEST_VCF), "GRCh38")
        return store

    def test_search_finds_variants(self, loaded_collection):
        samples = loaded_collection.list_samples()
        sample_table = samples[0][0]
        # Search a wide region around chr1 to find variants
        results = loaded_collection.search_region(sample_table, "chr1", 1, 1000000)
        assert isinstance(results, list)

    def test_search_empty_region(self, loaded_collection):
        samples = loaded_collection.list_samples()
        sample_table = samples[0][0]
        # Search in a region far from any variants in our small test VCF
        results = loaded_collection.search_region(sample_table, "chrY", 1000000, 2000000)
        assert results == [] or isinstance(results, list)

    def test_search_result_fields(self, loaded_collection):
        samples = loaded_collection.list_samples()
        sample_table = samples[0][0]
        # Wide search to find at least one result (within chr1 bounds)
        results = loaded_collection.search_region(sample_table, "chr1", 1, 248956422)
        if results:
            r = results[0]
            assert "uvid" in r
            assert "allele1" in r
            assert "allele2" in r
            assert "phased" in r
            assert "dp" in r
            assert "gq" in r
            assert "qual" in r
            assert "filter" in r
            assert "multiallelic" in r

    def test_search_invalid_chr(self, loaded_collection):
        samples = loaded_collection.list_samples()
        sample_table = samples[0][0]
        with pytest.raises(ValueError, match="Invalid chromosome"):
            loaded_collection.search_region(sample_table, "chr99", 1, 100)


# ──────────────────────────────────────────────────────
# CLI tests (using subprocess via typer testing)
# ──────────────────────────────────────────────────────

runner = CliRunner()


class TestCLIEncode:
    """Test the encode CLI command."""

    def test_encode_basic(self):
        result = runner.invoke(app, ["encode", "chr1", "100", "A", "G"])
        assert result.exit_code == 0
        assert "UVID:" in result.output
        assert "Integer:" in result.output

    def test_encode_with_assembly(self):
        result = runner.invoke(app, ["encode", "chr1", "100", "A", "G", "--assembly", "GRCh37"])
        assert result.exit_code == 0


class TestCLIDecode:
    """Test the decode CLI command."""

    def test_decode_basic(self):
        # First encode to get a valid hex
        uvid = UVID.encode("chr1", 100, "A", "G")
        hex_str = uvid.to_hex()

        result = runner.invoke(app, ["decode", hex_str])
        assert result.exit_code == 0
        assert "Chromosome:" in result.output
        assert "Position:" in result.output
        assert "REF:" in result.output
        assert "ALT:" in result.output
        assert "Assembly:" in result.output


class TestCLIAddAndQuery:
    """Test add, info, samples, and search CLI commands."""

    @pytest.fixture
    def collection_path(self, tmp_path):
        return str(tmp_path / "test.uvid")

    def test_add_vcf(self, collection_path):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        result = runner.invoke(app, ["add", collection_path, str(TEST_VCF)])
        assert result.exit_code == 0
        assert "Done." in result.output

    def test_info(self, collection_path):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        runner.invoke(app, ["add", collection_path, str(TEST_VCF)])
        result = runner.invoke(app, ["info", collection_path])
        assert result.exit_code == 0
        assert "Collection:" in result.output
        assert "Source files:" in result.output

    def test_samples(self, collection_path):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        runner.invoke(app, ["add", collection_path, str(TEST_VCF)])
        result = runner.invoke(app, ["samples", collection_path])
        assert result.exit_code == 0
        assert "Table" in result.output or "No samples" in result.output

    def test_search(self, collection_path):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        runner.invoke(app, ["add", collection_path, str(TEST_VCF)])

        # Get the sample table name, then close the collection before the CLI
        # reopens it (Windows DuckDB holds an exclusive file lock).
        store = Collection(collection_path)
        samples = store.list_samples()
        if not samples:
            pytest.skip("No samples in test VCF")
        sample_table = samples[0][0]
        del store

        result = runner.invoke(
            app,
            [
                "search",
                collection_path,
                "--sample",
                sample_table,
                "--chr",
                "chr1",
                "--start",
                "1",
                "--end",
                "248956422",
            ],
        )
        assert result.exit_code == 0

    def test_add_nonexistent_vcf(self, collection_path):
        result = runner.invoke(app, ["add", collection_path, "/nonexistent.vcf"])
        assert result.exit_code != 0

    def test_info_nonexistent_collection(self):
        result = runner.invoke(app, ["info", "/nonexistent.uvid"])
        assert result.exit_code != 0

    def test_samples_nonexistent_collection(self):
        result = runner.invoke(app, ["samples", "/nonexistent.uvid"])
        assert result.exit_code != 0


# ──────────────────────────────────────────────────────
# VCF passthrough tests
# ──────────────────────────────────────────────────────

SAMPLE_VCF = Path(__file__).parent.parent.parent / "test_data" / "sample.vcf"

# A minimal VCF with assembly in the header (for auto-detection tests)
VCF_WITH_ASSEMBLY = """\
##fileformat=VCFv4.3
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50
chr1\t200\trs456\tC\tT\t45\tPASS\tDP=60
"""


class TestVcfPassthrough:
    """Test the Rust-backed vcf_passthrough function."""

    def test_basic_with_assembly_override(self, tmp_path):
        """Passthrough with explicit assembly override on sample.vcf."""
        out = tmp_path / "out.vcf"
        count = vcf_passthrough(SAMPLE_VCF, out, assembly="GRCh38")
        assert count == 5

        lines = out.read_text().splitlines()
        data_lines = [l for l in lines if not l.startswith("#")]
        assert len(data_lines) == 5

        # First record: chr1 100 A>G — ID should be a UVID hex
        fields = data_lines[0].split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "100"
        assert len(fields[2]) == 35  # UVID hex format
        assert fields[3] == "A"
        assert fields[4] == "G"

    def test_assembly_auto_detection(self, tmp_path):
        """Assembly auto-detected from ##reference= header line."""
        vcf_in = tmp_path / "input.vcf"
        vcf_in.write_text(VCF_WITH_ASSEMBLY)
        out = tmp_path / "out.vcf"

        count = vcf_passthrough(vcf_in, out)  # no assembly override
        assert count == 2

        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        fields = data_lines[0].split("\t")
        assert len(fields[2]) == 35  # UVID hex

    def test_assembly_not_detected_error(self, tmp_path):
        """Raises AssemblyNotDetectedError when header has no assembly and no override."""
        out = tmp_path / "out.vcf"
        with pytest.raises(AssemblyNotDetectedError):
            vcf_passthrough(SAMPLE_VCF, out)

    def test_assembly_not_detected_is_value_error(self, tmp_path):
        """AssemblyNotDetectedError is a subclass of ValueError."""
        out = tmp_path / "out.vcf"
        with pytest.raises(ValueError):
            vcf_passthrough(SAMPLE_VCF, out)

    def test_uuid_mode(self, tmp_path):
        """--uuid flag produces UUIDv5 format IDs."""
        out = tmp_path / "out.vcf"
        vcf_passthrough(SAMPLE_VCF, out, use_uuid=True, assembly="GRCh38")

        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        fields = data_lines[0].split("\t")
        # UUID format: 8-4-4-4-12 = 36 chars
        assert len(fields[2]) == 36
        assert fields[2].count("-") == 4

    def test_multiallelic_id(self, tmp_path):
        """Multi-allelic records get underscore-joined IDs."""
        out = tmp_path / "out.vcf"
        vcf_passthrough(SAMPLE_VCF, out, assembly="GRCh38")

        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        # Record 4 (0-indexed 3) is chr2 1000 G A,C — multi-allelic
        fields = data_lines[3].split("\t")
        assert fields[0] == "chr2"
        assert fields[1] == "1000"
        assert "_" in fields[2]
        parts = fields[2].split("_")
        assert len(parts) == 2
        assert all(len(p) == 35 for p in parts)

    def test_dot_alt_keeps_dot_id(self, tmp_path):
        """ALT='.' records keep ID as '.'."""
        vcf_in = tmp_path / "dot_alt.vcf"
        vcf_in.write_text(
            "##fileformat=VCFv4.3\n"
            "##reference=GRCh38\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\n"
            "chr1\t200\t.\tC\t.\t10\tPASS\tDP=10\n"
        )
        out = tmp_path / "out.vcf"
        vcf_passthrough(vcf_in, out)

        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        # First record: chr1 100 A>G — should have a UVID
        fields0 = data_lines[0].split("\t")
        assert len(fields0[2]) == 35

        # Second record: chr1 200 C . — ALT is ".", ID stays "."
        fields1 = data_lines[1].split("\t")
        assert fields1[4] == "."
        assert fields1[2] == "."

    def test_header_passthrough_unchanged(self, tmp_path):
        """All header lines pass through byte-for-byte."""
        out = tmp_path / "out.vcf"
        vcf_passthrough(SAMPLE_VCF, out, assembly="GRCh38")

        original_headers = [l for l in SAMPLE_VCF.read_text().splitlines() if l.startswith("#")]
        output_headers = [l for l in out.read_text().splitlines() if l.startswith("#")]
        assert original_headers == output_headers

    def test_preserves_all_columns(self, tmp_path):
        """All columns except ID pass through unchanged."""
        out = tmp_path / "out.vcf"
        vcf_passthrough(SAMPLE_VCF, out, assembly="GRCh38")

        original_data = [l for l in SAMPLE_VCF.read_text().splitlines() if not l.startswith("#")]
        output_data = [l for l in out.read_text().splitlines() if not l.startswith("#")]

        for orig, out_line in zip(original_data, output_data, strict=False):
            orig_fields = orig.split("\t")
            out_fields = out_line.split("\t")
            # Columns 0, 1, 3+ should be identical
            assert orig_fields[0] == out_fields[0]  # CHROM
            assert orig_fields[1] == out_fields[1]  # POS
            assert orig_fields[3] == out_fields[3]  # REF
            assert orig_fields[4] == out_fields[4]  # ALT
            assert orig_fields[5] == out_fields[5]  # QUAL
            assert orig_fields[6] == out_fields[6]  # FILTER
            assert orig_fields[7] == out_fields[7]  # INFO
            if len(orig_fields) > 8:
                assert orig_fields[8] == out_fields[8]  # FORMAT
                assert orig_fields[9] == out_fields[9]  # SAMPLE1

    def test_deterministic(self, tmp_path):
        """Two runs produce identical output."""
        out1 = tmp_path / "out1.vcf"
        out2 = tmp_path / "out2.vcf"
        vcf_passthrough(SAMPLE_VCF, out1, assembly="GRCh38")
        vcf_passthrough(SAMPLE_VCF, out2, assembly="GRCh38")
        assert out1.read_text() == out2.read_text()

    def test_bgzf_output(self, tmp_path):
        """Output ending in .vcf.gz is bgzf-compressed."""
        out = tmp_path / "out.vcf.gz"
        count = vcf_passthrough(SAMPLE_VCF, out, assembly="GRCh38")
        assert count == 5

        # Check gzip magic bytes
        data = out.read_bytes()
        assert data[0:2] == b"\x1f\x8b"

    def test_stdout_output(self, tmp_path, capsys):
        """output=None writes to stdout."""
        # This writes to actual stdout, which capsys won't capture for Rust.
        # Instead, just verify it doesn't error.
        count = vcf_passthrough(SAMPLE_VCF, None, assembly="GRCh38")
        assert count == 5

    def test_uvid_matches_direct_encode(self, tmp_path):
        """UVIDs in passthrough output match UVID.encode() for the same variant."""
        out = tmp_path / "out.vcf"
        vcf_passthrough(SAMPLE_VCF, out, assembly="GRCh38")

        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        # First record: chr1 100 A G
        fields = data_lines[0].split("\t")
        passthrough_uvid = fields[2]

        direct_uvid = UVID.encode("chr1", 100, "A", "G", "GRCh38")
        assert passthrough_uvid == direct_uvid.to_hex()

    def test_uuid_matches_direct_encode(self, tmp_path):
        """UUIDs in passthrough output match UVID.encode().uuid5() for the same variant."""
        out = tmp_path / "out.vcf"
        vcf_passthrough(SAMPLE_VCF, out, use_uuid=True, assembly="GRCh38")

        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        fields = data_lines[0].split("\t")
        passthrough_uuid = fields[2]

        direct_uvid = UVID.encode("chr1", 100, "A", "G", "GRCh38")
        assert passthrough_uuid == str(direct_uvid.uuid5())


class TestCLIVcf:
    """Test the 'uvid vcf' CLI command."""

    def test_vcf_basic(self, tmp_path):
        out = tmp_path / "out.vcf"
        result = runner.invoke(app, ["vcf", str(SAMPLE_VCF), str(out), "-a", "GRCh38"])
        assert result.exit_code == 0
        assert "Processed 5 records" in result.output or "Processed 5 records" in (
            result.stderr or ""
        )

        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        assert len(data_lines) == 5

    def test_vcf_uuid_flag(self, tmp_path):
        out = tmp_path / "out.vcf"
        result = runner.invoke(app, ["vcf", str(SAMPLE_VCF), str(out), "-a", "GRCh38", "--uuid"])
        assert result.exit_code == 0

        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        fields = data_lines[0].split("\t")
        assert len(fields[2]) == 36  # UUID format

    def test_vcf_no_assembly_error(self, tmp_path):
        """Without assembly flag on a VCF with no assembly info, should error."""
        out = tmp_path / "out.vcf"
        result = runner.invoke(app, ["vcf", str(SAMPLE_VCF), str(out)])
        assert result.exit_code != 0

    def test_vcf_auto_detect_assembly(self, tmp_path):
        """VCF with ##reference=GRCh38 works without -a flag."""
        vcf_in = tmp_path / "input.vcf"
        vcf_in.write_text(VCF_WITH_ASSEMBLY)
        out = tmp_path / "out.vcf"

        result = runner.invoke(app, ["vcf", str(vcf_in), str(out)])
        assert result.exit_code == 0

    def test_vcf_nonexistent_input(self, tmp_path):
        out = tmp_path / "out.vcf"
        result = runner.invoke(app, ["vcf", "/nonexistent.vcf", str(out)])
        assert result.exit_code != 0

    def test_vcf_bgzf_output(self, tmp_path):
        out = tmp_path / "out.vcf.gz"
        result = runner.invoke(app, ["vcf", str(SAMPLE_VCF), str(out), "-a", "GRCh38"])
        assert result.exit_code == 0
        data = out.read_bytes()
        assert data[0:2] == b"\x1f\x8b"
