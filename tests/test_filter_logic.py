from genome_pipeline.filter_variants import apply_objective_filters, consequence_matches, passes_quality
from genome_pipeline.schemas import VariantRecord


def test_objective_filter_retains_rare_matching_variant():
    variant = VariantRecord(
        chrom="17",
        pos=43044295,
        ref="A",
        alt="G",
        qual="99",
        filter="PASS",
        gene="BRCA1",
        consequence="missense_variant",
        impact="MODERATE",
        gnomad_af=0.001,
    )
    objective = {
        "genes": ["BRCA1"],
        "consequence_filters": ["missense"],
        "max_allele_frequency": 0.01,
    }

    retained = apply_objective_filters([variant], objective, {"min_qual": 20, "require_pass_filter": True})

    assert retained == [variant]
    assert "gene matches objective" in retained[0].retention_reasons[0]


def test_objective_filter_excludes_common_variant():
    variant = VariantRecord(
        chrom="1",
        pos=1,
        ref="A",
        alt="G",
        qual="99",
        filter="PASS",
        gene="BRCA1",
        consequence="missense_variant",
        gnomad_af=0.2,
    )
    objective = {"genes": ["BRCA1"], "consequence_filters": ["missense"], "max_allele_frequency": 0.01}

    assert apply_objective_filters([variant], objective, {"min_qual": 20}) == []


def test_quality_and_consequence_helpers():
    variant = VariantRecord("1", 10, "A", "T", qual="30", filter="PASS", consequence="stop_gained&splice_region_variant")

    assert passes_quality(variant)
    assert consequence_matches(variant, ["splice"])
