"""Shared data structures for structured variant records."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


def normalize_chrom(chrom: str) -> str:
    """Return a chromosome label normalized for cross-table joins."""
    value = chrom.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value.upper()


def variant_key(chrom: str, pos: int | str, ref: str, alt: str) -> str:
    return f"{normalize_chrom(chrom)}:{int(pos)}:{ref.upper()}:{alt.upper()}"


@dataclass(slots=True)
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    id: str | None = None
    qual: str | None = None
    filter: str | None = None
    info: dict[str, Any] | None = None
    format_keys: tuple[str, ...] = ()
    sample_values: dict[str, str] | None = None
    genotype: str | None = None
    zygosity: str | None = None
    gene: str | None = None
    transcript: str | None = None
    consequence: str | None = None
    impact: str | None = None
    protein_change: str | None = None
    rsid: str | None = None
    clinvar_significance: str | None = None
    clinvar_trait: str | None = None
    clinvar_review_status: str | None = None
    clinvar_accession: str | None = None
    gnomad_af: float | None = None
    gnomad_source: str | None = None
    retention_reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    source: str | None = None

    @property
    def key(self) -> str:
        return variant_key(self.chrom, self.pos, self.ref, self.alt)

    def to_dict(self) -> dict[str, Any]:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "id": self.id,
            "rsid": self.rsid,
            "qual": self.qual,
            "filter": self.filter,
            "genotype": self.genotype,
            "zygosity": self.zygosity,
            "gene": self.gene,
            "transcript": self.transcript,
            "consequence": self.consequence,
            "impact": self.impact,
            "protein_change": self.protein_change,
            "clinvar_significance": self.clinvar_significance,
            "clinvar_trait": self.clinvar_trait,
            "clinvar_review_status": self.clinvar_review_status,
            "clinvar_accession": self.clinvar_accession,
            "gnomad_af": self.gnomad_af,
            "gnomad_source": self.gnomad_source,
            "retention_reasons": "; ".join(self.retention_reasons),
            "warnings": "; ".join(self.warnings),
            "source": self.source,
        }
