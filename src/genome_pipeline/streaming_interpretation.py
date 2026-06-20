"""Streaming interpretation outputs for large annotated gVCF-style WGS inputs."""

from __future__ import annotations

import heapq
from collections import Counter
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from .annotate import iter_vcf_records
from .clinvar import annotate_variant_clinvar, load_clinvar_table
from .filter_variants import objective_retention_reasons, passes_quality, variant_priority_key
from .gnomad import annotate_variant_gnomad, load_gnomad_table
from .io_utils import write_csv_rows, write_json, write_variant_outputs
from .schemas import VariantRecord
from .summarize import summarize_reduced_output
from .wgs_overview import generate_wgs_overview


def _normalized_label(value: str | None) -> str:
    if not value:
        return ""
    return value.strip().replace(" ", "_").lower()


def _is_pathogenic_label(value: str | None) -> bool:
    normalized = _normalized_label(value)
    return "pathogenic" in normalized and "conflicting" not in normalized


def _copy_variant(variant: VariantRecord, *, retention_reasons: list[str] | None = None) -> VariantRecord:
    return replace(
        variant,
        retention_reasons=list(retention_reasons or variant.retention_reasons),
        warnings=list(variant.warnings),
    )


def _push_top_variant(
    bucket: list[tuple[tuple[int, str, int], int, VariantRecord]],
    variant: VariantRecord,
    counter: int,
    limit: int,
) -> None:
    score = variant_priority_key(variant)
    entry = (score, counter, variant)
    if len(bucket) < limit:
        heapq.heappush(bucket, entry)
        return
    if (score, counter) > (bucket[0][0], bucket[0][1]):
        heapq.heapreplace(bucket, entry)


def _ordered_variants(bucket: list[tuple[tuple[int, str, int], int, VariantRecord]]) -> list[VariantRecord]:
    return [entry[2] for entry in sorted(bucket, key=lambda item: (item[0], item[1]), reverse=True)]


def _write_counter_csv(path: Path, key_name: str, counter: Counter[str]) -> None:
    items = sorted(counter.items(), key=lambda item: (-item[1], item[0]))
    write_csv_rows(path, [key_name, "count"], [{key_name: key, "count": count} for key, count in items])


def _write_markdown(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _count_consequences(counter: Counter[str], consequence: str | None) -> None:
    if not consequence:
        return
    for chunk in consequence.split(","):
        for term in chunk.split("&"):
            value = term.strip()
            if value:
                counter[value] += 1


@dataclass
class ObjectiveStreamState:
    label: str
    stem: str
    title: str
    output_dir: Path
    objective: dict[str, Any] | None
    limit: int
    retained_count: int = 0
    gene_counts: Counter[str] = field(default_factory=Counter)
    consequence_counts: Counter[str] = field(default_factory=Counter)
    impact_counts: Counter[str] = field(default_factory=Counter)
    clinvar_counts: Counter[str] = field(default_factory=Counter)
    bucket: list[tuple[tuple[int, str, int], int, VariantRecord]] = field(default_factory=list)

    def consider(self, variant: VariantRecord, filter_config: dict[str, Any], counter: int) -> None:
        if self.objective is None:
            return
        reasons = objective_retention_reasons(variant, self.objective, filter_config)
        if reasons is None:
            return

        self.retained_count += 1
        self.gene_counts[variant.gene or "unknown"] += 1
        if variant.impact:
            self.impact_counts[variant.impact.upper()] += 1
        if variant.clinvar_significance:
            self.clinvar_counts[variant.clinvar_significance] += 1
        _count_consequences(self.consequence_counts, variant.consequence)

        candidate = _copy_variant(variant, retention_reasons=reasons)
        _push_top_variant(self.bucket, candidate, counter, self.limit)

    @property
    def written_count(self) -> int:
        return len(self.bucket)


def _write_objective_outputs(
    state: ObjectiveStreamState,
    filter_config: dict[str, Any],
    logger,
) -> None:
    state.output_dir.mkdir(parents=True, exist_ok=True)
    if state.objective is None:
        _write_markdown(
            state.output_dir / "summary.md",
            [
                f"# {state.label}",
                "",
                "No disease-risk objective was requested for this run.",
            ],
        )
        write_json(
            state.output_dir / "summary.json",
            {
                "label": state.label,
                "objective_requested": False,
                "filter_config": filter_config,
            },
        )
        return

    ranked = _ordered_variants(state.bucket)
    write_variant_outputs(state.output_dir, state.stem, ranked, state.title)
    _write_counter_csv(state.output_dir / "gene_counts.csv", "gene", state.gene_counts)
    _write_counter_csv(state.output_dir / "consequence_counts.csv", "consequence", state.consequence_counts)
    _write_counter_csv(state.output_dir / "impact_counts.csv", "impact", state.impact_counts)
    _write_counter_csv(state.output_dir / "clinvar_significance_counts.csv", "clinvar_significance", state.clinvar_counts)

    summary_payload = {
        "label": state.label,
        "objective_requested": True,
        "objective_name": state.objective.get("name"),
        "workflow": state.objective.get("workflow"),
        "description": state.objective.get("description"),
        "configured_gene_count": len(state.objective.get("genes", [])),
        "consequence_filters": state.objective.get("consequence_filters", []),
        "max_allele_frequency": state.objective.get("max_allele_frequency"),
        "retained_variant_count": state.retained_count,
        "rows_written": len(ranked),
        "row_limit": state.limit,
        "truncated": state.retained_count > len(ranked),
    }
    write_json(state.output_dir / "summary.json", summary_payload)

    lines = [
        f"# {state.label}",
        "",
        "Research-only prioritization. Not medical advice or diagnosis.",
        "",
        state.objective.get("description", ""),
        "",
        f"Configured genes: {len(state.objective.get('genes', []))}",
        f"Retained variants matching the objective: {state.retained_count}",
        f"Rows written to `{state.stem}.csv`: {len(ranked)}",
        f"Row limit: {state.limit}",
        "",
        "Streaming mode was used to keep large WGS outputs compact and AI-friendly.",
    ]
    if state.retained_count > len(ranked):
        lines.extend(
            [
                "",
                "More variants matched than were written to the candidate table. Review the count tables and summary JSON for the full retained totals.",
            ]
        )
    _write_markdown(state.output_dir / "summary.md", lines)
    logger.info(
        "Wrote %s objective outputs with %s retained variants and %s rows exported",
        state.label,
        state.retained_count,
        len(ranked),
    )


def run_streaming_interpretation_pipeline(
    annotated_vcf: Path,
    normalized_vcf: Path,
    run_dir: Path,
    selected: dict[str, object],
    secondary_objective: dict[str, object],
    pgx_objective: dict[str, object],
    config: dict[str, object],
    logger,
) -> None:
    annotation_config = config.get("annotation", {})
    filter_config = config.get("filters", {})
    sample_index = int(annotation_config.get("sample_index", 0))
    max_records = annotation_config.get("max_records")
    preview_limit = int(annotation_config.get("preview_limit", 200))
    top_variant_limit = int(annotation_config.get("top_variant_limit", 100))
    candidate_output_limit = int(annotation_config.get("candidate_output_limit", 500))

    disease_objective = selected if selected.get("workflow") == "disease_risk" else None
    states = [
        ObjectiveStreamState(
            label="Disease Risk",
            stem="top_candidates",
            title="Disease Risk Top Candidates",
            output_dir=run_dir / "disease_risk",
            objective=disease_objective,
            limit=candidate_output_limit,
        ),
        ObjectiveStreamState(
            label="Secondary Findings",
            stem="secondary_findings",
            title="Secondary Findings Candidates",
            output_dir=run_dir / "secondary_findings",
            objective=secondary_objective,
            limit=candidate_output_limit,
        ),
        ObjectiveStreamState(
            label="Pharmacogenomics",
            stem="pgx_gene_candidates",
            title="Pharmacogenomics Gene Candidates",
            output_dir=run_dir / "pharmacogenomics",
            objective=pgx_objective,
            limit=candidate_output_limit,
        ),
    ]

    generate_wgs_overview(normalized_vcf, run_dir, config, logger)

    clinvar_table = load_clinvar_table(Path(config["resources"]["clinvar_table"])) if config.get("resources", {}).get("clinvar_table") else {}
    gnomad_table = load_gnomad_table(Path(config["resources"]["gnomad_table"])) if config.get("resources", {}).get("gnomad_table") else {}
    logger.info("Loaded %s local ClinVar rows for streaming interpretation", len(clinvar_table))
    logger.info("Loaded %s local gnomAD rows for streaming interpretation", len(gnomad_table))

    gene_counts: Counter[str] = Counter()
    consequence_counts: Counter[str] = Counter()
    impact_counts: Counter[str] = Counter()
    clinvar_counts: Counter[str] = Counter()
    pathogenic_preview: list[tuple[tuple[int, str, int], int, VariantRecord]] = []
    annotated_preview: list[VariantRecord] = []

    scanned = 0
    quality_passing = 0
    with_gene = 0
    with_consequence = 0
    with_impact = 0
    with_clinvar = 0
    with_gnomad = 0

    for variant in iter_vcf_records(
        annotated_vcf,
        sample_index=sample_index,
        max_records=max_records,
        parse_info_fields=True,
        skip_reference_blocks=True,
        skip_homozygous_reference=True,
    ):
        scanned += 1
        annotate_variant_clinvar(variant, clinvar_table)
        annotate_variant_gnomad(variant, gnomad_table)

        if variant.gene:
            with_gene += 1
        if variant.consequence:
            with_consequence += 1
        if variant.impact:
            with_impact += 1
        if variant.clinvar_significance:
            with_clinvar += 1
        if variant.gnomad_af is not None:
            with_gnomad += 1

        if not passes_quality(
            variant,
            min_qual=filter_config.get("min_qual"),
            require_pass_filter=filter_config.get("require_pass_filter", True),
        ):
            continue

        quality_passing += 1
        if len(annotated_preview) < preview_limit:
            annotated_preview.append(_copy_variant(variant))

        gene_counts[variant.gene or "unknown"] += 1
        if variant.impact:
            impact_counts[variant.impact.upper()] += 1
        if variant.clinvar_significance:
            clinvar_counts[variant.clinvar_significance] += 1
        _count_consequences(consequence_counts, variant.consequence)

        if _is_pathogenic_label(variant.clinvar_significance):
            _push_top_variant(pathogenic_preview, _copy_variant(variant), scanned, top_variant_limit)

        for state in states:
            state.consider(variant, filter_config, scanned)

    write_variant_outputs(
        run_dir / "annotated",
        "variants_with_refs_preview",
        annotated_preview,
        "Annotated Variants With Local References Preview",
    )
    write_variant_outputs(
        run_dir / "annotated",
        "clinvar_pathogenic_preview",
        _ordered_variants(pathogenic_preview),
        "ClinVar Pathogenic Or Likely Pathogenic Preview",
    )
    _write_counter_csv(run_dir / "annotated" / "gene_counts.csv", "gene", gene_counts)
    _write_counter_csv(run_dir / "annotated" / "consequence_counts.csv", "consequence", consequence_counts)
    _write_counter_csv(run_dir / "annotated" / "impact_counts.csv", "impact", impact_counts)
    _write_counter_csv(run_dir / "annotated" / "clinvar_significance_counts.csv", "clinvar_significance", clinvar_counts)

    annotation_summary = {
        "input_vcf": str(annotated_vcf),
        "normalized_vcf": str(normalized_vcf),
        "streaming_mode": True,
        "preview_limit": preview_limit,
        "top_variant_limit": top_variant_limit,
        "candidate_output_limit": candidate_output_limit,
        "max_records": max_records,
        "non_reference_annotated_alleles_scanned": scanned,
        "quality_passing_annotated_alleles": quality_passing,
        "with_gene_annotation": with_gene,
        "with_consequence_annotation": with_consequence,
        "with_impact_annotation": with_impact,
        "with_clinvar_annotation": with_clinvar,
        "with_gnomad_af": with_gnomad,
        "objective_retained_counts": {
            "disease_risk": states[0].retained_count,
            "secondary_findings": states[1].retained_count,
            "pharmacogenomics": states[2].retained_count,
        },
        "output_row_counts": {
            "annotated_preview": len(annotated_preview),
            "clinvar_pathogenic_preview": len(pathogenic_preview),
            "disease_risk": states[0].written_count,
            "secondary_findings": states[1].written_count,
            "pharmacogenomics": states[2].written_count,
        },
    }
    write_json(run_dir / "annotated" / "annotation_summary.json", annotation_summary)
    _write_markdown(
        run_dir / "annotated" / "summary.md",
        [
            "# Annotation Summary",
            "",
            "Research-only output. Not medical advice or diagnosis.",
            "",
            "Streaming mode was used because this is a large annotated gVCF-style WGS input.",
            "Full all-variant exports were intentionally skipped to keep outputs usable for people and downstream AI review.",
            "",
            f"Non-reference annotated alleles scanned: {scanned}",
            f"Quality-passing annotated alleles: {quality_passing}",
            f"Variants with gene annotations: {with_gene}",
            f"Variants with local ClinVar info: {with_clinvar}",
            f"Variants with gnomAD allele frequency: {with_gnomad}",
            "",
            f"Annotated preview rows written: {len(annotated_preview)}",
            f"Genome-wide ClinVar pathogenic preview rows written: {len(pathogenic_preview)}",
            "",
            "Use the objective directories for ranked candidate tables and the count CSVs for a compact genome-wide view.",
        ],
    )

    for state in states:
        _write_objective_outputs(state, filter_config, logger)

    summary_input = run_dir / "disease_risk" / "top_candidates.csv"
    if summary_input.exists():
        summarize_reduced_output(summary_input, run_dir / "summaries" / "disease_risk_summary.md")
