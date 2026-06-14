"""Local gnomAD frequency comparison interface."""

from __future__ import annotations

import csv
from pathlib import Path

from .schemas import VariantRecord, variant_key


def _as_float(value: object) -> float | None:
    if value in (None, "", "."):
        return None
    try:
        return float(str(value).split(",")[0])
    except ValueError:
        return None


def load_gnomad_table(path: Path | None) -> dict[str, dict[str, str]]:
    if not path or not path.exists():
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"chrom", "pos", "ref", "alt"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(f"gnomAD table {path} must include columns: {', '.join(sorted(required))}")
        return {
            variant_key(row["chrom"], row["pos"], row["ref"], row["alt"]): row
            for row in reader
            if row.get("chrom") and row.get("pos") and row.get("ref") and row.get("alt")
        }


def apply_existing_frequency_info(variant: VariantRecord) -> None:
    if variant.gnomad_af is not None:
        return
    info = variant.info or {}
    for key in ("gnomAD_AF", "GNOMAD_AF", "gnomADg_AF", "AF_POPMAX", "AF"):
        af = _as_float(info.get(key))
        if af is not None:
            variant.gnomad_af = af
            variant.gnomad_source = key
            return


def annotate_variant_gnomad(
    variant: VariantRecord,
    table: dict[str, dict[str, str]],
) -> bool:
    apply_existing_frequency_info(variant)
    row = table.get(variant.key)
    if not row:
        return False

    af = _as_float(row.get("af_popmax"))
    if af is None:
        af = _as_float(row.get("af"))
    if af is not None:
        variant.gnomad_af = af
        variant.gnomad_source = row.get("source") or "local_gnomad_table"
        return True
    return False


def compare_gnomad(variants: list[VariantRecord], table_path: Path | None, logger=None) -> list[VariantRecord]:
    table = load_gnomad_table(table_path)
    if logger:
        logger.info("Loaded %s local gnomAD rows from %s", len(table), table_path or "no table")

    matches = 0
    for variant in variants:
        if annotate_variant_gnomad(variant, table):
            matches += 1
    if logger:
        logger.info("Matched %s variants to local gnomAD data", matches)
    return variants


def is_rare(variant: VariantRecord, threshold: float | None, include_unknown: bool = True) -> bool:
    if threshold is None:
        return True
    if variant.gnomad_af is None:
        return include_unknown
    return variant.gnomad_af <= threshold
