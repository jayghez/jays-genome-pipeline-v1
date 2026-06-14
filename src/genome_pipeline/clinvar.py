"""Local ClinVar comparison interface."""

from __future__ import annotations

import csv
from pathlib import Path

from .schemas import VariantRecord, variant_key


def load_clinvar_table(path: Path | None) -> dict[str, dict[str, str]]:
    if not path or not path.exists():
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"chrom", "pos", "ref", "alt"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(f"ClinVar table {path} must include columns: {', '.join(sorted(required))}")
        return {
            variant_key(row["chrom"], row["pos"], row["ref"], row["alt"]): row
            for row in reader
            if row.get("chrom") and row.get("pos") and row.get("ref") and row.get("alt")
        }


def apply_existing_clinvar_info(variant: VariantRecord) -> None:
    info = variant.info or {}
    if not variant.clinvar_significance and info.get("CLNSIG"):
        variant.clinvar_significance = str(info["CLNSIG"])
    if not variant.clinvar_trait and info.get("CLNDN"):
        variant.clinvar_trait = str(info["CLNDN"]).replace("_", " ")
    if not variant.clinvar_review_status and info.get("CLNREVSTAT"):
        variant.clinvar_review_status = str(info["CLNREVSTAT"]).replace("_", " ")
    if not variant.clinvar_accession and info.get("ALLELEID"):
        variant.clinvar_accession = f"ALLELEID:{info['ALLELEID']}"


def annotate_variant_clinvar(
    variant: VariantRecord,
    table: dict[str, dict[str, str]],
) -> bool:
    apply_existing_clinvar_info(variant)
    row = table.get(variant.key)
    if not row:
        return False

    variant.clinvar_significance = row.get("clinvar_significance") or variant.clinvar_significance
    variant.clinvar_trait = row.get("disease_trait") or row.get("trait") or variant.clinvar_trait
    variant.clinvar_review_status = row.get("review_status") or variant.clinvar_review_status
    variant.clinvar_accession = row.get("accession") or variant.clinvar_accession
    return True


def compare_clinvar(variants: list[VariantRecord], table_path: Path | None, logger=None) -> list[VariantRecord]:
    table = load_clinvar_table(table_path)
    if logger:
        logger.info("Loaded %s local ClinVar rows from %s", len(table), table_path or "no table")

    matches = 0
    for variant in variants:
        if annotate_variant_clinvar(variant, table):
            matches += 1
    if logger:
        logger.info("Matched %s variants to local ClinVar data", matches)
    return variants
