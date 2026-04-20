"""Variant filtering and prioritization helpers."""

from __future__ import annotations

from typing import Any, Iterable

from .gnomad import is_rare
from .schemas import VariantRecord


CONSEQUENCE_GROUPS = {
    "protein_altering": {
        "missense_variant",
        "frameshift_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "inframe_insertion",
        "inframe_deletion",
        "protein_altering_variant",
    },
    "splice": {
        "splice_acceptor_variant",
        "splice_donor_variant",
        "splice_region_variant",
    },
    "missense": {"missense_variant"},
    "nonsense_frameshift": {"stop_gained", "frameshift_variant"},
}

PATHOGENIC_TERMS = {"pathogenic", "likely_pathogenic", "likely pathogenic"}


def passes_quality(variant: VariantRecord, min_qual: float | None = 20, require_pass_filter: bool = True) -> bool:
    if require_pass_filter and variant.filter not in (None, "PASS", "."):
        return False
    if variant.zygosity == "homozygous_ref":
        return False
    if min_qual is None or variant.qual in (None, "."):
        return True
    try:
        return float(variant.qual) >= min_qual
    except ValueError:
        return True


def expand_consequence_filters(filters: Iterable[str]) -> set[str]:
    terms: set[str] = set()
    for item in filters:
        if item in CONSEQUENCE_GROUPS:
            terms.update(CONSEQUENCE_GROUPS[item])
        else:
            terms.add(item)
    return terms


def consequence_matches(variant: VariantRecord, filters: Iterable[str]) -> bool:
    terms = expand_consequence_filters(filters)
    if not terms:
        return True
    if not variant.consequence:
        return False
    observed = {part.strip() for chunk in variant.consequence.split(",") for part in chunk.split("&")}
    return bool(observed & terms)


def gene_matches(variant: VariantRecord, genes: Iterable[str]) -> bool:
    wanted = {gene.upper() for gene in genes}
    if not wanted:
        return True
    return bool(variant.gene and variant.gene.upper() in wanted)


def apply_basic_filters(variants: Iterable[VariantRecord], filter_config: dict[str, Any]) -> list[VariantRecord]:
    return [
        variant
        for variant in variants
        if passes_quality(
            variant,
            min_qual=filter_config.get("min_qual"),
            require_pass_filter=filter_config.get("require_pass_filter", True),
        )
    ]


def apply_objective_filters(
    variants: Iterable[VariantRecord],
    objective: dict[str, Any],
    filter_config: dict[str, Any],
) -> list[VariantRecord]:
    retained: list[VariantRecord] = []
    genes = objective.get("genes", [])
    consequence_filters = objective.get("consequence_filters", [])
    threshold = objective.get("max_allele_frequency", filter_config.get("default_rare_af", 0.01))
    include_unknown_af = filter_config.get("include_unknown_af", True)

    for variant in variants:
        reasons: list[str] = []
        if not passes_quality(
            variant,
            min_qual=filter_config.get("min_qual"),
            require_pass_filter=filter_config.get("require_pass_filter", True),
        ):
            continue
        if genes:
            if not gene_matches(variant, genes):
                continue
            reasons.append(f"gene matches objective ({variant.gene})")
        if consequence_filters:
            if not consequence_matches(variant, consequence_filters):
                continue
            reasons.append(f"consequence matches objective ({variant.consequence})")
        if not is_rare(variant, threshold, include_unknown=include_unknown_af):
            continue
        if threshold is not None:
            if variant.gnomad_af is None:
                reasons.append("allele frequency unavailable; retained for review")
            else:
                reasons.append(f"allele frequency <= {threshold:g}")
        if variant.clinvar_significance:
            reasons.append(f"ClinVar: {variant.clinvar_significance}")
        variant.retention_reasons = reasons
        retained.append(variant)
    return retained


def rank_variants(variants: Iterable[VariantRecord]) -> list[VariantRecord]:
    def score(variant: VariantRecord) -> tuple[int, str, int]:
        value = 0
        sig = (variant.clinvar_significance or "").lower().replace(" ", "_")
        if "pathogenic" in sig and "conflicting" not in sig:
            value += 100 if "likely" not in sig else 80
        if variant.impact:
            value += {"HIGH": 40, "MODERATE": 25, "LOW": 5}.get(variant.impact.upper(), 0)
        if variant.gnomad_af is not None and variant.gnomad_af <= 0.01:
            value += 15
        if variant.gene:
            value += 5
        return (value, variant.chrom, -variant.pos)

    return sorted(variants, key=score, reverse=True)
