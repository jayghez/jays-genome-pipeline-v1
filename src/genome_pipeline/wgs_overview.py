"""Streaming WGS overview outputs for large raw gVCF-style inputs."""

from __future__ import annotations

import heapq
from collections import Counter
from pathlib import Path
from typing import Any

from .annotate import iter_vcf_records
from .filter_variants import passes_quality
from .io_utils import open_text_maybe_gzip, write_csv_rows, write_json, write_variant_outputs
from .schemas import VariantRecord, normalize_chrom


QUALITY_BINS = ["missing", "<20", "20-49", "50-99", "100+"]
HEADER_KEYS = [
    "fileDate",
    "source",
    "dataSourceType",
    "dataAnalysisProvider",
    "reference",
    "referenceInfo",
    "PipelineVersion",
]


def _qual_as_float(qual: str | None) -> float | None:
    if qual in (None, "", "."):
        return None
    try:
        return float(qual)
    except ValueError:
        return None


def classify_variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) == len(alt):
        return "mnv"
    if len(ref) != len(alt):
        return "indel"
    return "complex"


def quality_bin(qual: str | None) -> str:
    value = _qual_as_float(qual)
    if value is None:
        return "missing"
    if value < 20:
        return "<20"
    if value < 50:
        return "20-49"
    if value < 100:
        return "50-99"
    return "100+"


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    value = normalize_chrom(chrom)
    if value.isdigit():
        return (0, f"{int(value):02d}")
    if value == "X":
        return (1, value)
    if value == "Y":
        return (2, value)
    if value in {"M", "MT"}:
        return (3, value)
    return (4, value)


def is_primary_chrom(chrom: str) -> bool:
    value = normalize_chrom(chrom)
    return value in {str(index) for index in range(1, 23)} | {"X", "Y", "M", "MT"}


def _read_header_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with open_text_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                break
            if not line.startswith("##"):
                continue
            payload = line[2:].rstrip("\n")
            if "=" not in payload:
                continue
            key, value = payload.split("=", 1)
            if key in HEADER_KEYS:
                metadata[key] = value
    return metadata


def _push_top_variant(
    bucket: list[tuple[tuple[float, float, int], int, VariantRecord]],
    variant: VariantRecord,
    score: tuple[float, float, int],
    counter: int,
    limit: int,
) -> None:
    entry = (score, counter, variant)
    if len(bucket) < limit:
        heapq.heappush(bucket, entry)
        return
    if score > bucket[0][0]:
        heapq.heapreplace(bucket, entry)


def _ordered_variants(bucket: list[tuple[tuple[float, float, int], int, VariantRecord]]) -> list[VariantRecord]:
    return [entry[2] for entry in sorted(bucket, key=lambda item: (item[0], item[1]), reverse=True)]


def _write_counter_csv(path: Path, key_name: str, counter: Counter[str], sort_key=None) -> None:
    items = list(counter.items())
    if sort_key is None:
        items.sort(key=lambda item: item[0])
    else:
        items.sort(key=sort_key)
    write_csv_rows(path, [key_name, "count"], [{key_name: key, "count": count} for key, count in items])


def _write_markdown(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def generate_wgs_overview(
    normalized_vcf: Path,
    run_dir: Path,
    config: dict[str, Any],
    logger,
) -> dict[str, Any]:
    annotation_config = config.get("annotation", {})
    sample_index = int(annotation_config.get("sample_index", 0))
    max_records = annotation_config.get("max_records")
    preview_limit = int(annotation_config.get("preview_limit", 200))
    top_limit = int(annotation_config.get("top_variant_limit", 100))
    filter_config = config.get("filters", {})

    overview_dir = run_dir / "wgs_overview"
    header_metadata = _read_header_metadata(normalized_vcf)
    logger.info("Generating WGS overview outputs from %s", normalized_vcf)

    chromosome_counts: Counter[str] = Counter()
    variant_type_counts: Counter[str] = Counter()
    zygosity_counts: Counter[str] = Counter()
    filter_counts: Counter[str] = Counter()
    quality_counts: Counter[str] = Counter()

    variant_preview: list[VariantRecord] = []
    basic_preview: list[VariantRecord] = []
    top_homozygous_alt: list[tuple[tuple[float, float, int], int, VariantRecord]] = []
    top_no_rsid: list[tuple[tuple[float, float, int], int, VariantRecord]] = []
    top_long_indels: list[tuple[tuple[float, float, int], int, VariantRecord]] = []
    top_quality: list[tuple[tuple[float, float, int], int, VariantRecord]] = []

    non_reference = 0
    basic_quality = 0
    with_rsid = 0
    without_rsid = 0
    primary_assembly = 0
    counter = 0

    for variant in iter_vcf_records(
        normalized_vcf,
        sample_index=sample_index,
        max_records=max_records,
        parse_info_fields=False,
        skip_reference_blocks=True,
        skip_homozygous_reference=True,
    ):
        counter += 1
        non_reference += 1
        chromosome_counts[normalize_chrom(variant.chrom)] += 1
        variant_kind = classify_variant_type(variant.ref, variant.alt)
        variant_type_counts[variant_kind] += 1
        zygosity_counts[variant.zygosity or "unknown"] += 1
        filter_counts[variant.filter or "PASS"] += 1
        quality_counts[quality_bin(variant.qual)] += 1

        if variant.rsid:
            with_rsid += 1
        else:
            without_rsid += 1
        is_primary = is_primary_chrom(variant.chrom)
        if is_primary:
            primary_assembly += 1

        if len(variant_preview) < preview_limit:
            variant_preview.append(variant)

        qual_value = _qual_as_float(variant.qual) or -1.0
        indel_size = abs(len(variant.ref) - len(variant.alt))

        if passes_quality(
            variant,
            min_qual=filter_config.get("min_qual"),
            require_pass_filter=filter_config.get("require_pass_filter", True),
        ):
            basic_quality += 1
            if len(basic_preview) < preview_limit:
                basic_preview.append(variant)

            if is_primary:
                _push_top_variant(top_quality, variant, (qual_value, indel_size, -variant.pos), counter, top_limit)
                if variant.zygosity == "homozygous_alt":
                    _push_top_variant(top_homozygous_alt, variant, (qual_value, indel_size, -variant.pos), counter, top_limit)
                if not variant.rsid:
                    _push_top_variant(top_no_rsid, variant, (qual_value, indel_size, -variant.pos), counter, top_limit)
                if variant_kind == "indel":
                    _push_top_variant(top_long_indels, variant, (float(indel_size), qual_value, -variant.pos), counter, top_limit)

    pass_rate = round((basic_quality / non_reference), 4) if non_reference else 0.0
    rsid_rate = round((with_rsid / non_reference), 4) if non_reference else 0.0

    stats = {
        "input_metadata": header_metadata,
        "normalized_vcf": str(normalized_vcf),
        "has_structured_annotations": False,
        "has_gvcf_blocks": True,
        "annotation_backend_required": True,
        "non_reference_variant_alleles": non_reference,
        "basic_quality_variant_alleles": basic_quality,
        "basic_quality_pass_rate": pass_rate,
        "primary_assembly_variant_alleles": primary_assembly,
        "non_primary_assembly_variant_alleles": non_reference - primary_assembly,
        "with_rsid_count": with_rsid,
        "without_rsid_count": without_rsid,
        "with_rsid_rate": rsid_rate,
        "preview_limit": preview_limit,
        "top_variant_limit": top_limit,
        "max_records": max_records,
        "variant_type_counts": dict(variant_type_counts),
        "zygosity_counts": dict(zygosity_counts),
        "filter_counts": dict(filter_counts),
        "quality_bin_counts": {label: quality_counts.get(label, 0) for label in QUALITY_BINS},
        "chromosome_counts": dict(sorted(chromosome_counts.items(), key=lambda item: _chrom_sort_key(item[0]))),
    }

    write_json(overview_dir / "summary.json", stats)
    write_variant_outputs(run_dir / "annotated", "variant_preview", variant_preview, "Non-Reference Variant Preview")
    write_variant_outputs(run_dir / "filtered", "basic_quality_preview", basic_preview, "Basic Quality Variant Preview")
    write_variant_outputs(overview_dir, "top_pass_variants", _ordered_variants(top_quality), "Top PASS Variants By Quality")
    write_variant_outputs(
        overview_dir,
        "top_pass_homozygous_alt",
        _ordered_variants(top_homozygous_alt),
        "Top PASS Homozygous Alternate Variants",
    )
    write_variant_outputs(
        overview_dir,
        "top_pass_no_rsid",
        _ordered_variants(top_no_rsid),
        "Top PASS Variants Without rsID",
    )
    write_variant_outputs(
        overview_dir,
        "top_pass_long_indels",
        _ordered_variants(top_long_indels),
        "Top PASS Indels By Size",
    )

    _write_counter_csv(
        overview_dir / "chromosome_counts.csv",
        "chrom",
        chromosome_counts,
        sort_key=lambda item: _chrom_sort_key(item[0]),
    )
    _write_counter_csv(overview_dir / "variant_type_counts.csv", "variant_type", variant_type_counts)
    _write_counter_csv(overview_dir / "zygosity_counts.csv", "zygosity", zygosity_counts)
    _write_counter_csv(overview_dir / "filter_counts.csv", "filter", filter_counts, sort_key=lambda item: (-item[1], item[0]))
    _write_counter_csv(
        overview_dir / "quality_bin_counts.csv",
        "quality_bin",
        quality_counts,
        sort_key=lambda item: QUALITY_BINS.index(item[0]),
    )

    source = header_metadata.get("source", "unknown source")
    data_type = header_metadata.get("dataSourceType", "unknown data type")
    reference = header_metadata.get("referenceInfo") or header_metadata.get("reference", "unknown reference")
    provider = header_metadata.get("dataAnalysisProvider", "unknown provider")
    pipeline_version = header_metadata.get("PipelineVersion", "unknown")
    top_chromosomes = ", ".join(
        f"{chrom}: {count:,}"
        for chrom, count in sorted(chromosome_counts.items(), key=lambda item: item[1], reverse=True)[:5]
    )
    variant_breakdown = ", ".join(
        f"{kind}: {variant_type_counts.get(kind, 0):,}" for kind in ("snv", "indel", "mnv", "complex")
    )
    zygosity_breakdown = ", ".join(f"{kind}: {count:,}" for kind, count in sorted(zygosity_counts.items()))

    _write_markdown(
        overview_dir / "summary.md",
        [
            "# WGS Overview",
            "",
            "Research-only summary. Not medical advice or diagnosis.",
            "",
            f"Source: {source}",
            f"Data type: {data_type}",
            f"Provider: {provider}",
            f"Reference: {reference}",
            f"Pipeline version: {pipeline_version}",
            "",
            f"Non-reference variant alleles scanned: {non_reference:,}",
            f"Basic quality-pass alleles: {basic_quality:,} ({pass_rate:.1%})",
            f"Primary-assembly alleles: {primary_assembly:,}",
            f"Non-primary or decoy-contig alleles: {non_reference - primary_assembly:,}",
            f"Variants with rsID: {with_rsid:,} ({rsid_rate:.1%})",
            f"Variants without rsID: {without_rsid:,}",
            "",
            f"Variant types: {variant_breakdown}",
            f"Zygosity: {zygosity_breakdown}",
            f"Top chromosomes by variant count: {top_chromosomes}",
            "",
            "Useful reduced files for AI or downstream review:",
            "- `summary.json` for machine-readable counts and metadata.",
            "- `chromosome_counts.csv`, `variant_type_counts.csv`, and `zygosity_counts.csv` for compact distribution summaries.",
            "- `top_pass_variants.csv` for representative high-quality PASS variants on primary chromosomes.",
            "- `top_pass_homozygous_alt.csv` for strong homozygous alternate calls on primary chromosomes.",
            "- `top_pass_no_rsid.csv` for likely uncatalogued primary-chromosome variants.",
            "- `top_pass_long_indels.csv` for larger primary-chromosome indel events.",
            "",
            "Limitation: this input does not contain `ANN` or `CSQ` effect annotations, so these outputs describe the WGS callset but do not assign genes, consequences, or clinical interpretation.",
        ],
    )

    _write_markdown(
        overview_dir / "ai_context.md",
        [
            "# AI Context",
            "",
            "Use this file together with `summary.json` and the reduced CSV outputs.",
            "",
            f"This run comes from {source} on {reference}. It is a raw WGS gVCF-derived callset with {non_reference:,} non-reference alleles and {basic_quality:,} basic quality-pass alleles.",
            "",
            "What an AI can do well with these files:",
            "- Describe the scale and shape of the callset.",
            "- Compare chromosomes, variant classes, and zygosity distributions.",
            "- Review representative high-quality variants, homozygous alternate calls, and uncatalogued variants.",
            "- Help plan the next annotation or filtering step.",
            "",
            "What an AI should not do from these files alone:",
            "- Infer pathogenicity or disease risk from raw coordinates or rsIDs.",
            "- Claim gene-level effects without a real annotation backend such as VEP or SnpEff.",
            "- Treat these outputs as clinical findings.",
        ],
    )

    return stats
