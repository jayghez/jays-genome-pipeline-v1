"""Input validation, run-directory management, and output writers."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

from .schemas import VariantRecord


RUN_SUBDIRS = [
    "normalized",
    "annotated",
    "filtered",
    "wgs_overview",
    "disease_risk",
    "secondary_findings",
    "pharmacogenomics",
    "summaries",
    "logs",
]


@dataclass(frozen=True)
class VCFInput:
    path: Path
    compressed: bool
    index_path: Path | None


def validate_vcf_path(path: str | Path) -> VCFInput:
    vcf_path = Path(path).expanduser()
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF path does not exist: {vcf_path}")
    if not vcf_path.is_file():
        raise ValueError(f"VCF path is not a file: {vcf_path}")

    name = vcf_path.name.lower()
    compressed = name.endswith(".vcf.gz") or name.endswith(".vcf.bgz")
    uncompressed = name.endswith(".vcf")
    if not compressed and not uncompressed:
        raise ValueError("Input must end with .vcf, .vcf.gz, or .vcf.bgz")

    index_path = None
    if compressed:
        with vcf_path.open("rb") as handle:
            magic = handle.read(2)
        if magic != b"\x1f\x8b":
            raise ValueError(f"File has compressed suffix but is not gzip data: {vcf_path}")
        for candidate in (vcf_path.with_suffix(vcf_path.suffix + ".tbi"), vcf_path.with_suffix(vcf_path.suffix + ".csi")):
            if candidate.exists():
                index_path = candidate
                break

    return VCFInput(path=vcf_path, compressed=compressed, index_path=index_path)


def open_text_maybe_gzip(path: Path):
    if path.name.lower().endswith((".vcf.gz", ".vcf.bgz", ".gz")):
        return gzip.open(path, "rt")
    return path.open("rt")


def verify_vcf_header(path: Path) -> None:
    saw_fileformat = False
    saw_columns = False
    with open_text_maybe_gzip(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##fileformat=VCF"):
                saw_fileformat = True
            if line.startswith("#CHROM"):
                saw_columns = True
                break
            if line_number > 5000:
                break
    if not saw_fileformat or not saw_columns:
        raise ValueError(
            f"{path} does not look like a valid VCF. Expected ##fileformat and #CHROM header lines."
        )


def create_run_dir(base_dir: str | Path = "outputs") -> Path:
    base = Path(base_dir)
    base.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = base / f"run_{stamp}"
    suffix = 1
    while run_dir.exists():
        run_dir = base / f"run_{stamp}_{suffix}"
        suffix += 1
    for subdir in RUN_SUBDIRS:
        (run_dir / subdir).mkdir(parents=True, exist_ok=True)
    return run_dir


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv_rows(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown_table(path: Path, variants: list[VariantRecord], title: str, limit: int = 100) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [variant.to_dict() for variant in variants[:limit]]
    columns = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene",
        "consequence",
        "impact",
        "zygosity",
        "clinvar_significance",
        "gnomad_af",
        "retention_reasons",
    ]
    lines = [
        f"# {title}",
        "",
        "Research-only output. Not medical advice or diagnosis.",
        "",
        f"Variant rows: {len(variants)}",
        "",
    ]
    if not rows:
        lines.append("No variants retained for this output.")
    else:
        lines.append("| " + " | ".join(columns) + " |")
        lines.append("| " + " | ".join(["---"] * len(columns)) + " |")
        for row in rows:
            values = [str(row.get(column, "") if row.get(column, "") is not None else "") for column in columns]
            lines.append("| " + " | ".join(value.replace("|", "/") for value in values) + " |")
        if len(variants) > limit:
            lines.append("")
            lines.append(f"Showing first {limit} rows only. See CSV/JSON for full output.")
    path.write_text("\n".join(lines) + "\n")


def write_variant_outputs(output_dir: Path, stem: str, variants: Iterable[VariantRecord], title: str) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    variant_list = list(variants)
    rows = [variant.to_dict() for variant in variant_list]
    csv_path = output_dir / f"{stem}.csv"
    json_path = output_dir / f"{stem}.json"
    md_path = output_dir / f"{stem}.md"

    fieldnames = list(VariantRecord("0", 0, "N", "N").to_dict().keys())
    write_csv_rows(csv_path, fieldnames, rows)

    write_json(json_path, rows)
    write_markdown_table(md_path, variant_list, title)
    return {"csv": csv_path, "json": json_path, "markdown": md_path}


def read_gene_list(path: str | Path) -> list[str]:
    gene_path = Path(path)
    genes: list[str] = []
    for line in gene_path.read_text().splitlines():
        value = line.strip()
        if value and not value.startswith("#"):
            genes.append(value.upper())
    return genes
