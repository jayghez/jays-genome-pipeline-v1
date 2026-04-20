"""Command-line interface for the local-first genome pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path

from .annotate import annotate_vcf
from .clinvar import compare_clinvar
from .config import DEFAULT_OBJECTIVES_CONFIG, DEFAULT_PIPELINE_CONFIG, get_objective, load_objectives, load_pipeline_config, resource_path
from .filter_variants import apply_basic_filters
from .gnomad import compare_gnomad
from .io_utils import create_run_dir, validate_vcf_path, write_variant_outputs
from .logging_utils import configure_logging
from .normalize import normalize_vcf
from .pharmacogenomics import run_pharmacogenomics_pass
from .phenotype import run_disease_risk_pass
from .secondary_findings import run_secondary_findings_pass
from .summarize import summarize_reduced_output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genome-pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run the local VCF workflow")
    run_parser.add_argument("--vcf", required=True, help="Input .vcf, .vcf.gz, or .vcf.bgz")
    run_parser.add_argument("--objective", required=True, help="Objective name, for example hereditary_cancer or pgx")
    run_parser.add_argument("--config", default=str(DEFAULT_PIPELINE_CONFIG), help="Pipeline YAML config")
    run_parser.add_argument("--objectives", default=str(DEFAULT_OBJECTIVES_CONFIG), help="Objectives YAML config")
    run_parser.add_argument("--out-dir", default="outputs", help="Base output directory")
    run_parser.add_argument("--verbose", action="store_true", help="Enable debug logging")

    summarize_parser = subparsers.add_parser("summarize", help="Summarize a reduced output file")
    summarize_parser.add_argument("--input", required=True, help="Reduced CSV/JSON/Markdown input")
    summarize_parser.add_argument("--output", help="Markdown output path")
    summarize_parser.add_argument("--ai", action="store_true", help="Request optional AI summarization; disabled in v1")
    return parser


def run_pipeline(args: argparse.Namespace) -> Path:
    validate_vcf_path(args.vcf)
    config = load_pipeline_config(args.config)
    objectives = load_objectives(args.objectives)
    selected = get_objective(objectives, args.objective)
    secondary_objective = get_objective(objectives, "secondary_findings")
    pgx_objective = get_objective(objectives, "pgx")

    run_dir = create_run_dir(args.out_dir)
    logger = configure_logging(run_dir / "logs", verbose=args.verbose)
    logger.info("Created run directory: %s", run_dir)
    logger.info("Selected objective: %s", selected.get("name"))

    normalized_vcf = normalize_vcf(args.vcf, run_dir / "normalized", config, logger)
    variants = annotate_vcf(normalized_vcf, run_dir / "annotated", config, logger)
    compare_clinvar(variants, resource_path(config, "clinvar_table"), logger)
    compare_gnomad(variants, resource_path(config, "gnomad_table"), logger)
    write_variant_outputs(run_dir / "annotated", "variants_with_refs", variants, "Annotated Variants With Local References")

    filter_config = config.get("filters", {})
    basic_variants = apply_basic_filters(variants, filter_config)
    write_variant_outputs(run_dir / "filtered", "basic_quality", basic_variants, "Basic Quality Filtered Variants")

    disease_objective = selected if selected.get("workflow") == "disease_risk" else None
    run_disease_risk_pass(basic_variants, disease_objective, run_dir / "disease_risk", filter_config, logger)
    run_secondary_findings_pass(basic_variants, secondary_objective, run_dir / "secondary_findings", filter_config, logger)
    run_pharmacogenomics_pass(basic_variants, pgx_objective, run_dir / "pharmacogenomics", filter_config, logger)

    summary_input = run_dir / "disease_risk" / "top_candidates.csv"
    if summary_input.exists():
        summarize_reduced_output(summary_input, run_dir / "summaries" / "disease_risk_summary.md")
    logger.info("Pipeline complete: %s", run_dir)
    return run_dir


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "run":
        run_dir = run_pipeline(args)
        print(f"Run complete: {run_dir}")
        return 0
    if args.command == "summarize":
        output = summarize_reduced_output(args.input, args.output, use_ai=args.ai)
        print(f"Summary written: {output}")
        return 0
    parser.error("Unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
