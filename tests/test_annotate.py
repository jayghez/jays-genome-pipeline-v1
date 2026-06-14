from pathlib import Path

from genome_pipeline.annotate import inspect_annotation_headers, iter_vcf_records


def test_inspect_annotation_headers_detects_ann() -> None:
    vcf = Path("tests/fixtures/annotated_header.vcf")
    inspection = inspect_annotation_headers(vcf)

    assert inspection.has_structured_annotations is True
    assert inspection.has_gvcf_blocks is False
    assert "ANN" in inspection.annotation_fields


def test_iter_vcf_records_skips_reference_blocks_and_hom_ref(tmp_path: Path) -> None:
    vcf = tmp_path / "sample.g.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "1\t10\t.\tA\t.\t50\tPASS\tEND=20\tGT\t0/0\n"
        "1\t21\trs1\tC\tT\t50\tPASS\t.\tGT\t0/1\n"
        "1\t22\trs2\tG\tA\t50\tPASS\t.\tGT\t0/0\n"
    )

    variants = list(
        iter_vcf_records(
            vcf,
            parse_info_fields=False,
            skip_reference_blocks=True,
            skip_homozygous_reference=True,
        )
    )

    assert len(variants) == 1
    assert variants[0].pos == 21
    assert variants[0].alt == "T"
