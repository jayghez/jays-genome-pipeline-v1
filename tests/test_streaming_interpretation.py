import csv
import json
import logging
from pathlib import Path

from genome_pipeline.io_utils import create_run_dir
from genome_pipeline.streaming_interpretation import run_streaming_interpretation_pipeline


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_streaming_interpretation_writes_compact_candidate_outputs(tmp_path: Path) -> None:
    annotated_vcf = tmp_path / "annotated.g.vcf"
    annotated_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n"
        "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p'\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "1\t10\t.\tA\t.\t50\tPASS\tEND=20\tGT\t0/0\n"
        "1\t21\trs1\tC\tT,<NON_REF>\t100\tPASS\tANN=T|missense_variant|MODERATE|BRCA1|GENE1|transcript|NM_BRCA1|protein_coding|1/10|c.100C>T|p.Ala34Val\tGT\t0/1\n"
        "1\t22\trs2\tG\tA,<NON_REF>\t80\tPASS\tANN=A|stop_gained|HIGH|TP53|GENE2|transcript|NM_TP53|protein_coding|1/10|c.200G>A|p.Trp67*\tGT\t1/1\n"
        "1\t23\trs3\tT\tC,<NON_REF>\t90\tPASS\tANN=C|synonymous_variant|LOW|CYP2C19|GENE3|transcript|NM_CYP2C19|protein_coding|1/10|c.300T>C|p.=\tGT\t0/1\n"
        "1\t24\trs4\tA\tG,<NON_REF>\t10\tLowQual\tANN=G|missense_variant|MODERATE|BRCA1|GENE1|transcript|NM_BRCA1|protein_coding|1/10|c.101A>G|p.Lys34Arg\tGT\t0/1\n"
    )

    clinvar = tmp_path / "clinvar.csv"
    clinvar.write_text(
        "chrom,pos,ref,alt,clinvar_significance,disease_trait,review_status,accession\n"
        "1,22,G,A,Pathogenic,Example syndrome,reviewed by expert panel,VCV0001\n"
    )
    gnomad = tmp_path / "gnomad.csv"
    gnomad.write_text(
        "chrom,pos,ref,alt,af,af_popmax,source\n"
        "1,21,C,T,0.001,,test\n"
        "1,22,G,A,0.0,,test\n"
        "1,23,T,C,0.2,,test\n"
    )

    config = {
        "annotation": {
            "sample_index": 0,
            "max_records": None,
            "preview_limit": 10,
            "top_variant_limit": 10,
            "candidate_output_limit": 10,
        },
        "filters": {
            "min_qual": 20,
            "require_pass_filter": True,
            "include_unknown_af": True,
            "default_rare_af": 0.01,
        },
        "resources": {
            "clinvar_table": str(clinvar),
            "gnomad_table": str(gnomad),
        },
    }
    selected = {
        "name": "hereditary_cancer",
        "workflow": "disease_risk",
        "description": "Cancer-focused review",
        "genes": ["BRCA1", "TP53"],
        "consequence_filters": ["missense", "nonsense_frameshift"],
        "max_allele_frequency": 0.01,
    }
    secondary_objective = {
        "name": "secondary_findings",
        "workflow": "secondary_findings",
        "description": "Actionable gene review",
        "genes": ["TP53"],
        "consequence_filters": ["nonsense_frameshift"],
        "max_allele_frequency": 0.01,
    }
    pgx_objective = {
        "name": "pharmacogenomics",
        "workflow": "pharmacogenomics",
        "description": "PGx starter review",
        "genes": ["CYP2C19"],
        "consequence_filters": [],
        "max_allele_frequency": None,
    }

    run_dir = create_run_dir(tmp_path / "outputs")

    run_streaming_interpretation_pipeline(
        annotated_vcf,
        annotated_vcf,
        run_dir,
        selected,
        secondary_objective,
        pgx_objective,
        config,
        logging.getLogger("test"),
    )

    summary = json.loads((run_dir / "annotated" / "annotation_summary.json").read_text())
    disease_rows = _read_csv_rows(run_dir / "disease_risk" / "top_candidates.csv")
    secondary_rows = _read_csv_rows(run_dir / "secondary_findings" / "secondary_findings.csv")
    pgx_rows = _read_csv_rows(run_dir / "pharmacogenomics" / "pgx_gene_candidates.csv")
    pathogenic_rows = _read_csv_rows(run_dir / "annotated" / "clinvar_pathogenic_preview.csv")

    assert summary["non_reference_annotated_alleles_scanned"] == 4
    assert summary["quality_passing_annotated_alleles"] == 3
    assert summary["objective_retained_counts"]["disease_risk"] == 2
    assert summary["objective_retained_counts"]["secondary_findings"] == 1
    assert summary["objective_retained_counts"]["pharmacogenomics"] == 1
    assert {row["gene"] for row in disease_rows} == {"BRCA1", "TP53"}
    assert {row["gene"] for row in secondary_rows} == {"TP53"}
    assert {row["gene"] for row in pgx_rows} == {"CYP2C19"}
    assert {row["gene"] for row in pathogenic_rows} == {"TP53"}
    assert (run_dir / "summaries" / "disease_risk_summary.md").exists()
