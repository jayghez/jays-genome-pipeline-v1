"""Disease-risk or objective-focused interpretation pass."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .filter_variants import apply_objective_filters, rank_variants
from .io_utils import write_variant_outputs
from .schemas import VariantRecord


def run_disease_risk_pass(
    variants: list[VariantRecord],
    objective: dict[str, Any] | None,
    output_dir: Path,
    filter_config: dict[str, Any],
    logger,
) -> list[VariantRecord]:
    output_dir.mkdir(parents=True, exist_ok=True)
    if objective is None:
        (output_dir / "summary.md").write_text(
            "# Disease Risk\n\nNo disease-risk objective was requested for this run.\n"
        )
        return []

    logger.info("Running disease-risk objective: %s", objective.get("name"))
    retained = apply_objective_filters(variants, objective, filter_config)
    ranked = rank_variants(retained)
    write_variant_outputs(output_dir, "top_candidates", ranked, "Disease Risk Top Candidates")
    _write_summary(output_dir / "summary.md", objective, ranked)
    return ranked


def _write_summary(path: Path, objective: dict[str, Any], variants: list[VariantRecord]) -> None:
    lines = [
        f"# Disease Risk: {objective.get('name')}",
        "",
        "Research-only prioritization. Not medical advice or diagnosis.",
        "",
        objective.get("description", ""),
        "",
        f"Retained variants: {len(variants)}",
        "",
        "Known annotations and heuristic prioritization are kept separate. A retained row is a review candidate, not a clinical finding.",
    ]
    if not variants:
        lines.extend(
            [
                "",
                "No variants were retained. This can happen when the VCF lacks gene/consequence annotations, local ClinVar/gnomAD resources are empty, or no variants match the objective.",
            ]
        )
    path.write_text("\n".join(lines) + "\n")
