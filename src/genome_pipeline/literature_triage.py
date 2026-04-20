"""Prompt generation for optional literature triage after candidate reduction."""

from __future__ import annotations

from pathlib import Path


def build_literature_triage_prompt(reduced_input: str, template_path: str | Path = "configs/prompts/literature_triage.md") -> str:
    template = Path(template_path).read_text()
    return (
        f"{template}\n\n"
        "Reduced candidate input follows. Treat this as a triage aid only, not a diagnosis.\n\n"
        "```text\n"
        f"{reduced_input[:8000]}\n"
        "```\n"
    )


def write_literature_triage_prompt(input_path: str | Path, output_path: str | Path | None = None) -> Path:
    source = Path(input_path)
    if not source.exists():
        raise FileNotFoundError(f"Reduced input not found: {source}")
    output = Path(output_path) if output_path else source.with_suffix(source.suffix + ".literature_prompt.md")
    output.write_text(build_literature_triage_prompt(source.read_text()))
    return output
