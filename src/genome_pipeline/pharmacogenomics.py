"""Starter pharmacogenomics pass."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .filter_variants import apply_objective_filters, rank_variants
from .io_utils import write_variant_outputs
from .schemas import VariantRecord


def run_pharmacogenomics_pass(
    variants: list[VariantRecord],
    objective: dict[str, Any],
    output_dir: Path,
    filter_config: dict[str, Any],
    logger,
) -> list[VariantRecord]:
    logger.info("Running pharmacogenomics gene-focused pass")
    output_dir.mkdir(parents=True, exist_ok=True)
    retained = apply_objective_filters(variants, objective, filter_config)
    ranked = rank_variants(retained)
    write_variant_outputs(output_dir, "pgx_gene_candidates", ranked, "Pharmacogenomics Gene Candidates")
    _write_summary(output_dir / "summary.md", objective, ranked)
    return ranked


def _write_summary(path: Path, objective: dict[str, Any], variants: list[VariantRecord]) -> None:
    lines = [
        "# Pharmacogenomics",
        "",
        "Research-only output. Not medical advice or prescribing guidance.",
        "",
        "This v1 pass is gene-focused only. It does not call star alleles, diplotypes, metabolizer phenotypes, or medication recommendations.",
        "",
        f"Configured PGx genes: {len(objective.get('genes', []))}",
        f"Retained variants: {len(variants)}",
        "",
        "Future upgrade path: PharmCAT-style integration or another validated PGx caller.",
    ]
    path.write_text("\n".join(lines) + "\n")
