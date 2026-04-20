"""Configurable secondary-findings pass."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .filter_variants import apply_objective_filters, rank_variants
from .io_utils import write_variant_outputs
from .schemas import VariantRecord


def run_secondary_findings_pass(
    variants: list[VariantRecord],
    objective: dict[str, Any],
    output_dir: Path,
    filter_config: dict[str, Any],
    logger,
) -> list[VariantRecord]:
    logger.info("Running secondary-findings pass")
    output_dir.mkdir(parents=True, exist_ok=True)
    retained = apply_objective_filters(variants, objective, filter_config)
    ranked = rank_variants(retained)
    write_variant_outputs(output_dir, "secondary_findings", ranked, "Secondary Findings Candidates")
    _write_summary(output_dir / "summary.md", objective, ranked)
    return ranked


def _write_summary(path: Path, objective: dict[str, Any], variants: list[VariantRecord]) -> None:
    genes = objective.get("genes", [])
    lines = [
        "# Secondary Findings",
        "",
        "Research-only output. Not medical advice or diagnosis.",
        "",
        "This v1 pass uses a configurable starter actionable gene list. It is not a complete ACMG secondary findings implementation.",
        "",
        f"Configured genes: {len(genes)}",
        f"Retained variants: {len(variants)}",
        "",
        "Interpretation should be limited to whether a variant deserves further review. Do not infer actionability from this table alone.",
    ]
    path.write_text("\n".join(lines) + "\n")
