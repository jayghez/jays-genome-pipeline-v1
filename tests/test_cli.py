import argparse
import logging
from pathlib import Path
from types import SimpleNamespace

from genome_pipeline import cli


def _objectives() -> dict[str, dict[str, object]]:
    return {
        "hereditary_cancer": {
            "name": "hereditary_cancer",
            "workflow": "disease_risk",
            "genes": ["BRCA1", "TP53"],
        },
        "secondary_findings": {
            "name": "secondary_findings",
            "workflow": "secondary_findings",
            "genes": ["TP53"],
        },
        "pgx": {
            "name": "pgx",
            "workflow": "pharmacogenomics",
            "genes": ["CYP2C19"],
        },
    }


def test_run_pipeline_routes_structured_gvcf_to_streaming_mode(monkeypatch, tmp_path: Path) -> None:
    run_dir = tmp_path / "outputs" / "run-001"
    normalized_vcf = tmp_path / "normalized.vcf.gz"
    prepared_vcf = tmp_path / "prepared.vcf.gz"
    logger = logging.getLogger("test-cli-streaming")
    streaming_call: dict[str, object] = {}

    monkeypatch.setattr(cli, "validate_vcf_path", lambda path: None)
    monkeypatch.setattr(cli, "load_pipeline_config", lambda path: {"filters": {}, "annotation": {}})
    monkeypatch.setattr(cli, "load_objectives", lambda path: _objectives())
    monkeypatch.setattr(cli, "create_run_dir", lambda out_dir: run_dir)
    monkeypatch.setattr(cli, "configure_logging", lambda logs_dir, verbose=False: logger)
    monkeypatch.setattr(cli, "normalize_vcf", lambda vcf, output_dir, config, log: normalized_vcf)
    monkeypatch.setattr(cli, "annotate_effects", lambda path, output_dir, config, log: prepared_vcf)
    monkeypatch.setattr(
        cli,
        "inspect_annotation_headers",
        lambda path: SimpleNamespace(has_structured_annotations=True, has_gvcf_blocks=True),
    )

    def fake_streaming(
        annotated_vcf: Path,
        normalized_input: Path,
        pipeline_run_dir: Path,
        selected: dict[str, object],
        secondary: dict[str, object],
        pgx: dict[str, object],
        config: dict[str, object],
        log,
    ) -> None:
        streaming_call.update(
            {
                "annotated_vcf": annotated_vcf,
                "normalized_vcf": normalized_input,
                "run_dir": pipeline_run_dir,
                "selected": selected,
                "secondary": secondary,
                "pgx": pgx,
                "config": config,
                "logger": log,
            }
        )

    monkeypatch.setattr(cli, "run_streaming_interpretation_pipeline", fake_streaming)
    monkeypatch.setattr(cli, "annotate_vcf", lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("annotate_vcf should not run")))
    monkeypatch.setattr(
        cli,
        "_run_unstructured_input_pipeline",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("preview mode should not run")),
    )

    args = argparse.Namespace(
        vcf="input.vcf.gz",
        objective="hereditary_cancer",
        config="configs/pipeline.yaml",
        objectives="configs/objectives.yaml",
        out_dir=str(tmp_path / "outputs"),
        verbose=False,
    )

    result = cli.run_pipeline(args)

    assert result == run_dir
    assert streaming_call["annotated_vcf"] == prepared_vcf
    assert streaming_call["normalized_vcf"] == normalized_vcf
    assert streaming_call["run_dir"] == run_dir
    assert streaming_call["selected"]["name"] == "hereditary_cancer"
    assert streaming_call["secondary"]["name"] == "secondary_findings"
    assert streaming_call["pgx"]["name"] == "pgx"


def test_run_pipeline_routes_unstructured_input_to_preview_mode(monkeypatch, tmp_path: Path) -> None:
    run_dir = tmp_path / "outputs" / "run-001"
    normalized_vcf = tmp_path / "normalized.vcf.gz"
    prepared_vcf = tmp_path / "prepared.vcf.gz"
    logger = logging.getLogger("test-cli-preview")
    preview_call: dict[str, object] = {}

    monkeypatch.setattr(cli, "validate_vcf_path", lambda path: None)
    monkeypatch.setattr(cli, "load_pipeline_config", lambda path: {"filters": {}, "annotation": {}})
    monkeypatch.setattr(cli, "load_objectives", lambda path: _objectives())
    monkeypatch.setattr(cli, "create_run_dir", lambda out_dir: run_dir)
    monkeypatch.setattr(cli, "configure_logging", lambda logs_dir, verbose=False: logger)
    monkeypatch.setattr(cli, "normalize_vcf", lambda vcf, output_dir, config, log: normalized_vcf)
    monkeypatch.setattr(cli, "annotate_effects", lambda path, output_dir, config, log: prepared_vcf)
    monkeypatch.setattr(
        cli,
        "inspect_annotation_headers",
        lambda path: SimpleNamespace(has_structured_annotations=False, has_gvcf_blocks=True),
    )

    def fake_preview(
        normalized_input: Path,
        pipeline_run_dir: Path,
        selected: dict[str, object],
        secondary: dict[str, object],
        pgx: dict[str, object],
        config: dict[str, object],
        log,
    ) -> None:
        preview_call.update(
            {
                "normalized_vcf": normalized_input,
                "run_dir": pipeline_run_dir,
                "selected": selected,
                "secondary": secondary,
                "pgx": pgx,
                "config": config,
                "logger": log,
            }
        )

    monkeypatch.setattr(cli, "_run_unstructured_input_pipeline", fake_preview)
    monkeypatch.setattr(
        cli,
        "run_streaming_interpretation_pipeline",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("streaming mode should not run")),
    )
    monkeypatch.setattr(cli, "annotate_vcf", lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("annotate_vcf should not run")))

    args = argparse.Namespace(
        vcf="input.vcf.gz",
        objective="hereditary_cancer",
        config="configs/pipeline.yaml",
        objectives="configs/objectives.yaml",
        out_dir=str(tmp_path / "outputs"),
        verbose=False,
    )

    result = cli.run_pipeline(args)

    assert result == run_dir
    assert preview_call["normalized_vcf"] == normalized_vcf
    assert preview_call["run_dir"] == run_dir
    assert preview_call["selected"]["name"] == "hereditary_cancer"
    assert preview_call["secondary"]["name"] == "secondary_findings"
    assert preview_call["pgx"]["name"] == "pgx"
