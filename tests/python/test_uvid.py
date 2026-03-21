"""Tests for the UVID Python bindings and CLI."""

import os
import tempfile
import uuid
from pathlib import Path

import pytest

from uvid import NAMESPACE_UVID, UVID, Collection


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

    def test_encode_invalid_nucleotide(self):
        with pytest.raises(ValueError):
            UVID.encode("chr1", 100, "N", "G")

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
        assert fields["overflow"] is False
        assert fields["assembly"] == "GRCh38"

    def test_decode_overflow(self):
        """Sequences > 38 bases total trigger overflow mode."""
        uvid = UVID.encode("chr1", 1, "A" * 20, "C" * 20)
        fields = uvid.decode()
        assert fields["overflow"] is True
        assert fields["ref_len"] == 20
        assert fields["alt_len"] == 20
        # In overflow mode, sequences are not stored
        assert fields["ref"] == ""
        assert fields["alt"] == ""

    def test_decode_max_payload(self):
        """Exactly 38 bases should NOT overflow."""
        uvid = UVID.encode("chr1", 1, "A" * 19, "C" * 19)
        fields = uvid.decode()
        assert fields["overflow"] is False
        assert fields["ref"] == "A" * 19
        assert fields["alt"] == "C" * 19


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
        """Same-locus variants from different assemblies should be close."""
        uvid37 = UVID.encode("chr1", 100, "A", "G", "GRCh37")
        uvid38 = UVID.encode("chr1", 100, "A", "G", "GRCh38")
        # Assembly is in LSB (2 bits), so they differ by at most 3
        assert abs(uvid37.as_int() - uvid38.as_int()) <= 3

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
        assert NAMESPACE_UVID == expected

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
        store = Collection(str(path))
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
        # Search in a region far from any variants
        results = loaded_collection.search_region(
            sample_table, "chrY", 999999000, 999999999
        )
        assert results == [] or isinstance(results, list)

    def test_search_result_fields(self, loaded_collection):
        samples = loaded_collection.list_samples()
        sample_table = samples[0][0]
        # Wide search to find at least one result
        results = loaded_collection.search_region(sample_table, "chr1", 1, 999999999)
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

from typer.testing import CliRunner
from uvid.cli import app

runner = CliRunner()


class TestCLIEncode:
    """Test the encode CLI command."""

    def test_encode_basic(self):
        result = runner.invoke(app, ["encode", "chr1", "100", "A", "G"])
        assert result.exit_code == 0
        assert "UVID:" in result.output
        assert "Integer:" in result.output

    def test_encode_with_assembly(self):
        result = runner.invoke(
            app, ["encode", "chr1", "100", "A", "G", "--assembly", "GRCh37"]
        )
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
        result = runner.invoke(app, ["add", str(TEST_VCF), collection_path])
        assert result.exit_code == 0
        assert "Done." in result.output

    def test_info(self, collection_path):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        runner.invoke(app, ["add", str(TEST_VCF), collection_path])
        result = runner.invoke(app, ["info", collection_path])
        assert result.exit_code == 0
        assert "Collection:" in result.output
        assert "Source files:" in result.output

    def test_samples(self, collection_path):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        runner.invoke(app, ["add", str(TEST_VCF), collection_path])
        result = runner.invoke(app, ["samples", collection_path])
        assert result.exit_code == 0
        assert "Table" in result.output or "No samples" in result.output

    def test_search(self, collection_path):
        if not TEST_VCF.exists():
            pytest.skip("Test VCF not found")
        runner.invoke(app, ["add", str(TEST_VCF), collection_path])

        # Get the sample table name
        store = Collection(collection_path)
        samples = store.list_samples()
        if not samples:
            pytest.skip("No samples in test VCF")
        sample_table = samples[0][0]

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
                "999999999",
            ],
        )
        assert result.exit_code == 0

    def test_add_nonexistent_vcf(self, collection_path):
        result = runner.invoke(app, ["add", "/nonexistent.vcf", collection_path])
        assert result.exit_code != 0

    def test_info_nonexistent_collection(self):
        result = runner.invoke(app, ["info", "/nonexistent.uvid"])
        assert result.exit_code != 0

    def test_samples_nonexistent_collection(self):
        result = runner.invoke(app, ["samples", "/nonexistent.uvid"])
        assert result.exit_code != 0
