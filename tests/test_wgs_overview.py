from genome_pipeline.wgs_overview import classify_variant_type, is_primary_chrom, quality_bin


def test_classify_variant_type() -> None:
    assert classify_variant_type("A", "G") == "snv"
    assert classify_variant_type("AT", "GC") == "mnv"
    assert classify_variant_type("A", "AT") == "indel"


def test_quality_bin() -> None:
    assert quality_bin(None) == "missing"
    assert quality_bin("19.9") == "<20"
    assert quality_bin("20") == "20-49"
    assert quality_bin("50") == "50-99"
    assert quality_bin("100") == "100+"


def test_is_primary_chrom() -> None:
    assert is_primary_chrom("1") is True
    assert is_primary_chrom("chrX") is True
    assert is_primary_chrom("MT") is True
    assert is_primary_chrom("Un_KN707963v1_decoy") is False
