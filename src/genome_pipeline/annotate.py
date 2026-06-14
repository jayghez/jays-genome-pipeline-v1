"""Structured VCF parsing and lightweight annotation extraction."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

from .io_utils import open_text_maybe_gzip, verify_vcf_header, write_variant_outputs
from .schemas import VariantRecord


IMPACT_ORDER = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
DEFAULT_ANN_FIELDS = [
    "Allele",
    "Annotation",
    "Annotation_Impact",
    "Gene_Name",
    "Gene_ID",
    "Feature_Type",
    "Feature_ID",
    "Transcript_BioType",
    "Rank",
    "HGVS.c",
    "HGVS.p",
]


@dataclass(frozen=True)
class AnnotationInspection:
    annotation_fields: dict[str, list[str]]
    has_gvcf_blocks: bool = False

    @property
    def has_structured_annotations(self) -> bool:
        return bool(self.annotation_fields)


def parse_info(value: str) -> dict[str, Any]:
    if value in {"", "."}:
        return {}
    parsed: dict[str, Any] = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" not in item:
            parsed[item] = True
            continue
        key, raw = item.split("=", 1)
        parsed[key] = raw
    return parsed


def parse_annotation_headers(header_lines: list[str]) -> dict[str, list[str]]:
    fields: dict[str, list[str]] = {}
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            match = re.search(r"Format: ([^\">]+)", line)
            if match:
                fields["CSQ"] = [part.strip() for part in match.group(1).split("|")]
        if line.startswith("##INFO=<ID=ANN"):
            fields["ANN"] = DEFAULT_ANN_FIELDS
    return fields


def inspect_annotation_headers(path: Path) -> AnnotationInspection:
    verify_vcf_header(path)
    header_lines: list[str] = []
    has_gvcf_blocks = False
    with open_text_maybe_gzip(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##"):
                header_lines.append(line)
                if line.startswith("##INFO=<ID=END"):
                    has_gvcf_blocks = True
                if line.startswith("##ALT=<ID=NON_REF"):
                    has_gvcf_blocks = True
                continue
            if line.startswith("#CHROM"):
                break
    return AnnotationInspection(
        annotation_fields=parse_annotation_headers(header_lines),
        has_gvcf_blocks=has_gvcf_blocks,
    )


def parse_format(format_value: str, sample_value: str | None) -> tuple[list[str], dict[str, str]]:
    if not format_value or format_value == "." or not sample_value:
        return [], {}
    keys = format_value.split(":")
    values = sample_value.split(":")
    return keys, {key: values[index] if index < len(values) else "" for index, key in enumerate(keys)}


def extract_genotype(format_value: str, sample_value: str | None) -> str | None:
    if not format_value or format_value == "." or not sample_value:
        return None
    keys = format_value.split(":")
    try:
        gt_index = keys.index("GT")
    except ValueError:
        return None
    values = sample_value.split(":")
    if gt_index >= len(values):
        return None
    return values[gt_index]


def derive_zygosity(gt: str | None) -> str | None:
    if not gt or gt in {".", "./.", ".|."}:
        return None
    alleles = re.split(r"[/|]", gt)
    if all(allele == "0" for allele in alleles):
        return "homozygous_ref"
    called = [allele for allele in alleles if allele != "."]
    alt = [allele for allele in called if allele != "0"]
    if len(set(called)) == 1 and alt:
        return "homozygous_alt"
    if len(alt) == 1 and "0" in called:
        return "heterozygous"
    if len(set(alt)) > 1:
        return "multi_alt"
    if alt:
        return "non_ref"
    return None


def _info_first(info: dict[str, Any], keys: list[str]) -> str | None:
    for key in keys:
        value = info.get(key)
        if value not in (None, "", "."):
            return str(value).split(",")[0]
    return None


def _best_annotation(items: list[dict[str, str]]) -> dict[str, str]:
    if not items:
        return {}
    return max(items, key=lambda item: IMPACT_ORDER.get(item.get("IMPACT", item.get("Annotation_Impact", "")).upper(), 0))


def extract_structured_annotation(info: dict[str, Any], alt: str, annotation_fields: dict[str, list[str]]) -> dict[str, str | None]:
    result: dict[str, str | None] = {
        "gene": _info_first(info, ["SYMBOL", "GENE", "Gene", "Gene.refGene", "GN"]),
        "transcript": _info_first(info, ["Feature", "Transcript", "transcript", "NM"]),
        "consequence": _info_first(info, ["Consequence", "Func.refGene", "ExonicFunc.refGene"]),
        "impact": _info_first(info, ["IMPACT", "Annotation_Impact"]),
        "protein_change": _info_first(info, ["HGVSp", "AAChange.refGene", "Protein_Change"]),
    }

    if "CSQ" in info and "CSQ" in annotation_fields:
        annotations = []
        for raw in str(info["CSQ"]).split(","):
            values = raw.split("|")
            item = {field: values[index] if index < len(values) else "" for index, field in enumerate(annotation_fields["CSQ"])}
            if item.get("Allele") in {alt, ""}:
                annotations.append(item)
        best = _best_annotation(annotations)
        result.update(
            {
                "gene": best.get("SYMBOL") or best.get("Gene") or result["gene"],
                "transcript": best.get("Feature") or result["transcript"],
                "consequence": best.get("Consequence") or result["consequence"],
                "impact": best.get("IMPACT") or result["impact"],
                "protein_change": best.get("HGVSp") or result["protein_change"],
            }
        )

    if "ANN" in info:
        annotations = []
        for raw in str(info["ANN"]).split(","):
            values = raw.split("|")
            item = {field: values[index] if index < len(values) else "" for index, field in enumerate(DEFAULT_ANN_FIELDS)}
            if item.get("Allele") in {alt, ""}:
                annotations.append(item)
        best = _best_annotation(annotations)
        result.update(
            {
                "gene": best.get("Gene_Name") or result["gene"],
                "transcript": best.get("Feature_ID") or result["transcript"],
                "consequence": best.get("Annotation") or result["consequence"],
                "impact": best.get("Annotation_Impact") or result["impact"],
                "protein_change": best.get("HGVS.p") or result["protein_change"],
            }
        )

    return result


def iter_vcf_records(
    path: Path,
    sample_index: int = 0,
    max_records: int | None = None,
    *,
    parse_info_fields: bool = True,
    skip_reference_blocks: bool = False,
    skip_homozygous_reference: bool = False,
) -> Iterator[VariantRecord]:
    verify_vcf_header(path)
    header_lines: list[str] = []
    column_header: list[str] = []
    annotation_fields: dict[str, list[str]] = {}
    yielded = 0

    with open_text_maybe_gzip(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                column_header = line.lstrip("#").split("\t")
                annotation_fields = parse_annotation_headers(header_lines)
                continue
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, record_id, ref, alts, qual, filter_value, info_value = parts[:8]
            if skip_reference_blocks and alts == ".":
                continue
            format_value = parts[8] if len(parts) > 8 else ""
            sample_column_index = 9 + sample_index
            sample_value = parts[sample_column_index] if len(parts) > sample_column_index else None
            genotype = extract_genotype(format_value, sample_value)
            zygosity = derive_zygosity(genotype)
            if skip_homozygous_reference and zygosity == "homozygous_ref":
                continue
            info = parse_info(info_value) if parse_info_fields and info_value not in {"", "."} else None

            for alt in alts.split(","):
                if skip_reference_blocks and alt == ".":
                    continue
                structured = (
                    extract_structured_annotation(info or {}, alt, annotation_fields)
                    if parse_info_fields
                    else {"gene": None, "transcript": None, "consequence": None, "impact": None, "protein_change": None}
                )
                rsid = record_id if record_id.startswith("rs") else _info_first(info or {}, ["RS", "RSID", "dbSNP"])
                yield VariantRecord(
                    chrom=chrom,
                    pos=int(pos),
                    id=None if record_id == "." else record_id,
                    ref=ref,
                    alt=alt,
                    qual=None if qual == "." else qual,
                    filter=None if filter_value == "." else filter_value,
                    info=info,
                    genotype=genotype,
                    zygosity=zygosity,
                    gene=structured["gene"],
                    transcript=structured["transcript"],
                    consequence=structured["consequence"],
                    impact=structured["impact"],
                    protein_change=structured["protein_change"],
                    rsid=rsid,
                    source=str(path),
                )
                yielded += 1
                if max_records is not None and yielded >= max_records:
                    return
    if not column_header:
        raise ValueError(f"No #CHROM header found in {path}")


def parse_vcf_records(path: Path, sample_index: int = 0, max_records: int | None = None) -> list[VariantRecord]:
    return list(iter_vcf_records(path, sample_index=sample_index, max_records=max_records))


def annotate_vcf(path: Path, output_dir: Path, config: dict[str, Any], logger) -> list[VariantRecord]:
    logger.info("Parsing structured VCF annotations from %s", path)
    annotation_config = config.get("annotation", {})
    variants = parse_vcf_records(
        path,
        sample_index=int(annotation_config.get("sample_index", 0)),
        max_records=annotation_config.get("max_records"),
    )
    logger.info("Parsed %s variant alleles", len(variants))
    write_variant_outputs(output_dir, "variants", variants, "Annotated Variants")
    return variants
