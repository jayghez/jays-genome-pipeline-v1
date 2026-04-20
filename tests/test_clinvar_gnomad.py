from pathlib import Path

from genome_pipeline.clinvar import compare_clinvar
from genome_pipeline.gnomad import compare_gnomad, is_rare
from genome_pipeline.schemas import VariantRecord


def test_clinvar_local_table_match(tmp_path: Path):
    table = tmp_path / "clinvar.csv"
    table.write_text(
        "chrom,pos,ref,alt,clinvar_significance,disease_trait,review_status,accession\n"
        "chr1,100,A,G,Pathogenic,Example trait,reviewed by expert panel,VCV000001\n"
    )
    variant = VariantRecord("1", 100, "A", "G")

    compare_clinvar([variant], table)

    assert variant.clinvar_significance == "Pathogenic"
    assert variant.clinvar_accession == "VCV000001"


def test_gnomad_local_table_match_and_rarity(tmp_path: Path):
    table = tmp_path / "gnomad.csv"
    table.write_text("chrom,pos,ref,alt,af,af_popmax,source\n1,100,A,G,0.001,,test\n")
    variant = VariantRecord("chr1", 100, "A", "G")

    compare_gnomad([variant], table)

    assert variant.gnomad_af == 0.001
    assert is_rare(variant, 0.01)


def test_gnomad_zero_af_popmax_is_preserved(tmp_path: Path):
    table = tmp_path / "gnomad.csv"
    table.write_text("chrom,pos,ref,alt,af,af_popmax,source\n1,100,A,G,0.1,0.0,test\n")
    variant = VariantRecord("1", 100, "A", "G")

    compare_gnomad([variant], table)

    assert variant.gnomad_af == 0.0
