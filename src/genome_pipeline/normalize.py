"""VCF normalization through bcftools subprocess wrappers."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

from .io_utils import validate_vcf_path, verify_vcf_header, write_json


class MissingToolError(RuntimeError):
    pass


def require_tool(tool: str) -> str:
    path = shutil.which(tool)
    if path:
        return path
    raise MissingToolError(
        f"{tool} is required for normalization but was not found on PATH.\n"
        "Install with: brew install bcftools htslib\n"
        "Or with conda: conda install -c bioconda bcftools htslib"
    )


def run_command(args: list[str], logger) -> None:
    logger.info("Running: %s", " ".join(args))
    completed = subprocess.run(args, text=True, capture_output=True, check=False)
    if completed.stdout:
        logger.debug(completed.stdout.strip())
    if completed.stderr:
        logger.debug(completed.stderr.strip())
    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code "
            f"{completed.returncode}: {' '.join(args)}\n{completed.stderr.strip()}"
        )


def normalize_vcf(input_path: str | Path, output_dir: Path, config: dict[str, Any], logger) -> Path:
    """Validate, bgzip if needed, split multiallelics, and optionally left-align."""
    bcftools = require_tool("bcftools")
    vcf_input = validate_vcf_path(input_path)
    verify_vcf_header(vcf_input.path)
    output_dir.mkdir(parents=True, exist_ok=True)

    prepared_input = vcf_input.path
    if not vcf_input.compressed:
        prepared_input = output_dir / "input.bgz.vcf.gz"
        run_command([bcftools, "view", "-Oz", "-o", str(prepared_input), str(vcf_input.path)], logger)
        run_command([bcftools, "index", "--tbi", "--force", str(prepared_input)], logger)

    normalized_path = output_dir / "normalized.vcf.gz"
    norm_args = [bcftools, "norm"]
    reference_fasta = config.get("normalize", {}).get("reference_fasta")
    if reference_fasta:
        reference_path = Path(reference_fasta)
        if not reference_path.exists():
            raise FileNotFoundError(f"Configured reference FASTA does not exist: {reference_path}")
        norm_args.extend(["-f", str(reference_path)])
    else:
        logger.info("No reference FASTA configured; splitting multiallelics without left-alignment.")

    if config.get("normalize", {}).get("split_multiallelics", True):
        norm_args.append("-m-any")
    norm_args.extend(["-Oz", "-o", str(normalized_path), str(prepared_input)])
    run_command(norm_args, logger)
    run_command([bcftools, "index", "--tbi", "--force", str(normalized_path)], logger)
    verify_vcf_header(normalized_path)

    write_json(
        output_dir / "normalization_report.json",
        {
            "input": str(vcf_input.path),
            "input_compressed": vcf_input.compressed,
            "input_index": str(vcf_input.index_path) if vcf_input.index_path else None,
            "normalized_vcf": str(normalized_path),
            "reference_fasta": reference_fasta,
            "split_multiallelics": config.get("normalize", {}).get("split_multiallelics", True),
        },
    )
    return normalized_path
