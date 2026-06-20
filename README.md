# Local-First Genome Pipeline

This repository is a Python-first starter workflow for exploring a personal or research VCF locally. It is built around a real compressed GIAB HG001 / NA12878 GRCh38 sample VCF so the path through the code looks like the path you would use later with a Sequencing.com-style `.vcf` or `.vcf.gz`.

It is for research and personal exploration only. It is not diagnostic software, not medical advice, and not a substitute for a clinician or certified genetic counselor.

## What It Does

- Validates `.vcf` and `.vcf.gz` inputs.
- Creates timestamped run folders under `outputs/`.
- Normalizes variants with `bcftools`, including multiallelic splitting and optional left alignment if a reference FASTA is configured.
- Can add SnpEff consequence annotations to raw VCF or gVCF inputs when a local SnpEff install is available.
- Parses structured VCF `INFO` and `FORMAT` fields into CSV, JSON, and Markdown.
- Runs separate passes for disease-risk objectives, secondary findings, and pharmacogenomics.
- Supports local, pluggable ClinVar and gnomAD comparison tables.
- Reduces candidate sets before any interpretation or optional AI summarization.
- Keeps AI off by default and never sends raw VCF data to an AI step.
- Detects large gVCF-style inputs without gene/effect annotations and writes preview plus readiness summaries instead of stalling on full all-variant exports.
- Streams large annotated gVCF-style WGS inputs into compact candidate tables, count summaries, and AI-friendly previews instead of writing massive all-variant exports.

## What It Does Not Do

- It does not diagnose disease.
- It does not make clinical pharmacogenomic diplotype or prescribing calls.
- It does not perform VEP annotation.
- It does not bundle ClinVar, gnomAD, reference FASTA, PharmCAT, or HPO databases.
- It does not decide pathogenicity for variants of uncertain significance.

## Setup

Install Python 3.10+ and `bcftools`. On macOS:

```bash
brew install bcftools htslib
```

Or with conda:

```bash
conda install -c bioconda bcftools htslib
```

Then install the Python package:

```bash
make setup
source .venv/bin/activate
make test
make check-tools
```

If your WGS file is a raw Sequencing.com-style gVCF without `ANN` or `CSQ`, also install Java 21+ and a local SnpEff copy under `resources/snpEff/`, then set or keep the default `annotation` values in `configs/pipeline.yaml`.

## Run The GIAB-Derived Sample

```bash
genome-pipeline run \
  --vcf data/sample/HG001_chr20_subset.vcf.gz \
  --objective hereditary_cancer
```

Other configured objective examples:

```bash
genome-pipeline run --vcf data/sample/HG001_chr20_subset.vcf.gz --objective cardiovascular_risk
genome-pipeline run --vcf data/sample/HG001_chr20_subset.vcf.gz --objective secondary_findings
genome-pipeline run --vcf data/sample/HG001_chr20_subset.vcf.gz --objective pgx
```

Each run writes:

```text
outputs/run_YYYYMMDD_HHMMSS/
  normalized/
  annotated/
  filtered/
  wgs_overview/
  disease_risk/
  secondary_findings/
  pharmacogenomics/
  summaries/
  logs/
```

## Swap In A Sequencing.com VCF

Put the file somewhere outside git, for example:

```text
data/private/my_sequencing_com_export.vcf.gz
```

Then run:

```bash
genome-pipeline run \
  --vcf data/private/my_sequencing_com_export.vcf.gz \
  --objective hereditary_cancer
```

The code accepts `.vcf` and `.vcf.gz`. Uncompressed VCFs are converted to bgzipped VCFs by `bcftools` during normalization. If you have a GRCh38 reference FASTA, set `normalize.reference_fasta` in `configs/pipeline.yaml` to enable left alignment against that reference.

If the input is a raw gVCF-style file without `ANN` or `CSQ` consequence annotations and no local SnpEff backend is available, the pipeline finishes by writing:

- normalized outputs
- non-reference and basic-quality preview tables
- WGS overview tables for chromosome counts, variant-type counts, zygosity counts, and representative high-quality variants
- annotation-readiness summaries that explain why disease-risk and PGx candidate tables are empty

If local SnpEff is available, the pipeline can annotate that file automatically and then continue in streaming mode. In that mode it writes compact genome-wide count tables plus ranked candidate outputs such as:

- `annotated/summary.md`
- `annotated/annotation_summary.json`
- `annotated/variants_with_refs_preview.csv`
- `annotated/clinvar_pathogenic_preview.csv`
- `annotated/gene_counts.csv`
- `disease_risk/top_candidates.csv`
- `secondary_findings/secondary_findings.csv`
- `pharmacogenomics/pgx_gene_candidates.csv`

## Local ClinVar And gnomAD Tables

Version 1 uses local CSV interfaces so the architecture is ready for real local joins later. Configure paths in `configs/pipeline.yaml`:

```yaml
resources:
  clinvar_table: configs/clinvar_stub.csv
  gnomad_table: configs/gnomad_stub.csv
```

Expected ClinVar columns:

```text
chrom,pos,ref,alt,clinvar_significance,disease_trait,review_status,accession
```

Expected gnomAD columns:

```text
chrom,pos,ref,alt,af,af_popmax,source
```

## AI Summaries

AI is optional and off by default. The CLI summary command consumes reduced CSV, JSON, or Markdown outputs:

```bash
genome-pipeline summarize \
  --input outputs/run_YYYYMMDD_HHMMSS/disease_risk/top_candidates.csv
```

Prompt templates live under `configs/prompts/`. They explicitly warn against inventing pathogenicity, overinterpreting VUS, or mixing weak and strong evidence.

## Privacy

VCFs are sensitive personal genomic data. Keep private files out of git, avoid cloud sync folders, and inspect outputs before sharing. The `.gitignore` intentionally excludes common VCF and generated output patterns.

## Upgrade Path

Planned future upgrades are marked with TODOs in the code:

- VEP annotation backend.
- Real local ClinVar VCF/SQLite joins.
- Real local gnomAD frequency joins.
- HPO/phenotype scoring.
- PharmCAT-style pharmacogenomics.
- Literature retrieval against local indexes or controlled APIs.
