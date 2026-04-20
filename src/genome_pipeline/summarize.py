"""Reduced-output summarization helpers. AI is optional and off by default."""

from __future__ import annotations

import csv
import json
from pathlib import Path


DISCLAIMER = (
    "Research-only summary. Not medical advice or diagnosis. Do not overinterpret VUS, "
    "and do not infer pathogenicity that is not present in the reduced input."
)


def _read_reduced_input(path: Path) -> tuple[int, str]:
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text())
        return (len(data) if isinstance(data, list) else 1, path.read_text()[:4000])
    if path.suffix.lower() == ".csv":
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        preview = "\n".join(str(row) for row in rows[:20])
        return len(rows), preview
    text = path.read_text()
    return len([line for line in text.splitlines() if line.strip()]), text[:4000]


def summarize_reduced_output(input_path: str | Path, output_path: str | Path | None = None, use_ai: bool = False) -> Path:
    path = Path(input_path)
    if not path.exists():
        raise FileNotFoundError(f"Reduced input not found: {path}")
    row_count, preview = _read_reduced_input(path)
    output = Path(output_path) if output_path else path.with_suffix(path.suffix + ".summary.md")

    if use_ai:
        raise RuntimeError(
            "AI summarization is intentionally disabled in v1. Use configs/prompts/cautious_summary.md "
            "with reduced outputs only after explicitly enabling an AI integration."
        )

    lines = [
        "# Reduced Output Summary",
        "",
        DISCLAIMER,
        "",
        f"Input: `{path}`",
        f"Approximate row/item count: {row_count}",
        "",
        "This deterministic summary does not reinterpret pathogenicity. It only confirms that a reduced file was produced.",
        "",
        "## Preview",
        "",
        "```text",
        preview or "No content.",
        "```",
    ]
    output.write_text("\n".join(lines) + "\n")
    return output
