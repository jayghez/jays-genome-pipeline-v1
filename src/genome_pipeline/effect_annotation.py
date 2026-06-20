"""External consequence annotation backends for raw VCF and gVCF inputs."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

from .annotate import inspect_annotation_headers
from .io_utils import open_text_maybe_gzip, verify_vcf_header, write_json
from .normalize import MissingToolError, require_tool


MIN_SNPEFF_JAVA_MAJOR = 21
DEFAULT_SNPEFF_DIR = Path("resources/snpEff")
DEFAULT_SNPEFF_GENOME = "GRCh38.p13"
COMMON_JAVA21_PATHS = [
    Path("/opt/homebrew/opt/openjdk@21/bin/java"),
    Path("/usr/local/opt/openjdk@21/bin/java"),
]


def parse_java_major_version(version_text: str) -> int | None:
    match = re.search(r'version "([^"]+)"', version_text)
    if not match:
        return None
    raw = match.group(1)
    parts = raw.split(".")
    if not parts:
        return None
    try:
        major = int(parts[0])
    except ValueError:
        return None
    if major == 1 and len(parts) > 1:
        try:
            return int(parts[1])
        except ValueError:
            return None
    return major


def _java_version(java_cmd: str) -> int | None:
    completed = subprocess.run([java_cmd, "-version"], text=True, capture_output=True, check=False)
    version_text = (completed.stderr or completed.stdout).strip()
    return parse_java_major_version(version_text)


def _resolve_java_cmd(config: dict[str, Any]) -> str:
    annotation_config = config.get("annotation", {})
    configured = annotation_config.get("java_cmd")
    candidates: list[str] = []
    if configured:
        candidates.append(str(configured))
    java_on_path = shutil.which("java")
    if java_on_path:
        candidates.append(java_on_path)
    for candidate in COMMON_JAVA21_PATHS:
        if candidate.exists():
            candidates.append(str(candidate))

    checked: list[tuple[str, int | None]] = []
    seen: set[str] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        major = _java_version(candidate)
        checked.append((candidate, major))
        if major is not None and major >= MIN_SNPEFF_JAVA_MAJOR:
            return candidate

    details = ", ".join(f"{cmd} (Java {major or 'unknown'})" for cmd, major in checked) or "no java executables found"
    raise MissingToolError(
        "SnpEff requires Java 21 or newer. Checked: "
        f"{details}. Install a newer JDK or set annotation.java_cmd."
    )


def _resolve_snpeff_paths(config: dict[str, Any]) -> tuple[Path, Path, Path]:
    annotation_config = config.get("annotation", {})
    jar_value = annotation_config.get("snpeff_jar")
    if jar_value:
        jar_path = Path(str(jar_value))
        config_path = jar_path.parent / "snpEff.config"
        data_dir = jar_path.parent / "data"
        return jar_path, config_path, data_dir

    base_dir = Path(str(annotation_config.get("snpeff_dir", DEFAULT_SNPEFF_DIR)))
    jar_path = base_dir / "snpEff.jar"
    config_path = Path(str(annotation_config.get("snpeff_config", base_dir / "snpEff.config")))
    data_dir = Path(str(annotation_config.get("snpeff_data_dir", base_dir / "data")))
    return jar_path, config_path, data_dir


def _snpeff_genome_name(config: dict[str, Any]) -> str:
    annotation_config = config.get("annotation", {})
    return str(annotation_config.get("snpeff_genome", DEFAULT_SNPEFF_GENOME))


def _read_header_reference(path: Path) -> str:
    verify_vcf_header(path)
    lines: list[str] = []
    with open_text_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                break
            if line.startswith("##"):
                lines.append(line.rstrip("\n"))
    for line in lines:
        if line.startswith("##referenceInfo="):
            return line.split("=", 1)[1]
    for line in lines:
        if line.startswith("##reference="):
            return line.split("=", 1)[1]
    for line in lines:
        if "assembly=human_GRCh38" in line:
            return "human_GRCh38_no_alt_analysis_set"
    return "unknown"


def annotate_effects(path: Path, output_dir: Path, config: dict[str, Any], logger) -> Path:
    inspection = inspect_annotation_headers(path)
    if inspection.has_structured_annotations:
        return path

    annotation_config = config.get("annotation", {})
    backend = str(annotation_config.get("backend", "auto")).lower()
    if backend not in {"auto", "none", "snpeff"}:
        raise ValueError(f"Unsupported annotation backend: {backend}")
    if backend == "none":
        return path

    try:
        return _annotate_with_snpeff(path, output_dir, config, logger)
    except Exception as exc:
        if backend == "snpeff":
            raise
        logger.warning("External annotation skipped: %s", exc)
        return path


def _annotate_with_snpeff(path: Path, output_dir: Path, config: dict[str, Any], logger) -> Path:
    java_cmd = _resolve_java_cmd(config)
    snpeff_jar, snpeff_config, data_dir = _resolve_snpeff_paths(config)
    if not snpeff_jar.exists():
        raise MissingToolError(
            f"SnpEff jar not found at {snpeff_jar}. Install SnpEff under resources/ or set annotation.snpeff_jar."
        )
    if not snpeff_config.exists():
        raise MissingToolError(
            f"SnpEff config not found at {snpeff_config}. Install SnpEff under resources/ or set annotation.snpeff_config."
        )

    bcftools = require_tool("bcftools")
    genome = _snpeff_genome_name(config)
    memory_gb = int(config.get("annotation", {}).get("snpeff_memory_gb", 8))
    extra_args = [str(arg) for arg in config.get("annotation", {}).get("snpeff_extra_args", [])]
    local_db = data_dir / genome / "snpEffectPredictor.bin"

    output_dir.mkdir(parents=True, exist_ok=True)
    annotated_vcf = output_dir / "snpeff_annotated.vcf.gz"
    snpeff_log = output_dir / "snpeff.stderr.log"
    bcftools_log = output_dir / "bcftools.stderr.log"

    snpeff_args = [
        java_cmd,
        f"-Xmx{memory_gb}g",
        "-jar",
        str(snpeff_jar),
        "-c",
        str(snpeff_config.resolve()),
        "-dataDir",
        str(data_dir.resolve()),
    ]
    snpeff_args.extend(
        [
            *(["-nodownload"] if local_db.exists() else []),
            "-noStats",
            *extra_args,
            genome,
            str(path),
        ]
    )
    bcftools_args = [bcftools, "view", "-Oz", "-o", str(annotated_vcf)]

    logger.info("Running functional annotation: %s | %s", " ".join(snpeff_args), " ".join(bcftools_args))
    with snpeff_log.open("w") as snpeff_handle, bcftools_log.open("w") as bcftools_handle:
        snpeff = subprocess.Popen(snpeff_args, stdout=subprocess.PIPE, stderr=snpeff_handle)
        assert snpeff.stdout is not None
        bcf_view = subprocess.Popen(bcftools_args, stdin=snpeff.stdout, stderr=bcftools_handle)
        snpeff.stdout.close()
        bcf_code = bcf_view.wait()
        snpeff_code = snpeff.wait()

    if snpeff_code != 0 or bcf_code != 0:
        snippet = (snpeff_log.read_text() + "\n" + bcftools_log.read_text()).strip()
        raise RuntimeError(
            "SnpEff annotation failed.\n"
            f"SnpEff exit code: {snpeff_code}\n"
            f"bcftools exit code: {bcf_code}\n"
            f"{snippet[-4000:]}"
        )

    subprocess.run([bcftools, "index", "--tbi", "--force", str(annotated_vcf)], check=True, capture_output=True, text=True)
    verify_vcf_header(annotated_vcf)
    inspection = inspect_annotation_headers(annotated_vcf)
    if not inspection.has_structured_annotations:
        raise RuntimeError(f"SnpEff completed but {annotated_vcf} still has no ANN or CSQ annotations")

    report = {
        "input_vcf": str(path),
        "annotated_vcf": str(annotated_vcf),
        "backend": "snpeff",
        "snpeff_jar": str(snpeff_jar),
        "snpeff_config": str(snpeff_config),
        "snpeff_data_dir": str(data_dir),
        "snpeff_genome": genome,
        "local_database_present": local_db.exists(),
        "java_cmd": java_cmd,
        "java_min_major": MIN_SNPEFF_JAVA_MAJOR,
        "java_major": _java_version(java_cmd),
        "input_reference": _read_header_reference(path),
        "snpeff_stderr_log": str(snpeff_log),
        "bcftools_stderr_log": str(bcftools_log),
    }
    write_json(output_dir / "annotation_report.json", report)
    return annotated_vcf
