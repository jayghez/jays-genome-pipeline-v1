from pathlib import Path

import pytest

from genome_pipeline.io_utils import create_run_dir, validate_vcf_path, verify_vcf_header


def test_validate_vcf_path_accepts_uncompressed_vcf(tmp_path: Path):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    result = validate_vcf_path(vcf)
    verify_vcf_header(vcf)

    assert result.path == vcf
    assert result.compressed is False


def test_validate_vcf_path_rejects_missing(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        validate_vcf_path(tmp_path / "missing.vcf")


def test_create_run_dir_has_expected_subdirs(tmp_path: Path):
    run_dir = create_run_dir(tmp_path)

    assert (run_dir / "normalized").is_dir()
    assert (run_dir / "logs").is_dir()
    assert (run_dir / "wgs_overview").is_dir()
