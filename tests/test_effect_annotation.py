import logging
from pathlib import Path

from genome_pipeline.effect_annotation import annotate_effects, parse_java_major_version


def test_parse_java_major_version_handles_modern_and_legacy_formats() -> None:
    assert parse_java_major_version('openjdk version "21.0.3" 2024-04-16') == 21
    assert parse_java_major_version('java version "1.8.0_402"') == 8
    assert parse_java_major_version("not a java version string") is None


def test_annotate_effects_returns_original_when_ann_is_already_present(tmp_path: Path) -> None:
    vcf = tmp_path / "annotated.vcf"
    vcf.write_text(Path("tests/fixtures/annotated_header.vcf").read_text() + "1\t100\t.\tA\tG\t50\tPASS\tANN=G|missense_variant|MODERATE|BRCA1|GENE1|transcript|NM_1|protein_coding|1/2|c.1A>G|p.Lys1Arg\n")

    result = annotate_effects(vcf, tmp_path / "out", {"annotation": {"backend": "auto"}}, logging.getLogger("test"))

    assert result == vcf
