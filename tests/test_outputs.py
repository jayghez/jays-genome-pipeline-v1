from pathlib import Path

from genome_pipeline.io_utils import write_variant_outputs
from genome_pipeline.schemas import VariantRecord
from genome_pipeline.summarize import summarize_reduced_output


def test_write_variant_outputs_and_summary(tmp_path: Path):
    variant = VariantRecord("1", 100, "A", "G", gene="BRCA1", consequence="missense_variant")

    paths = write_variant_outputs(tmp_path, "candidates", [variant], "Candidates")
    summary = summarize_reduced_output(paths["csv"])

    assert paths["csv"].exists()
    assert paths["json"].exists()
    assert paths["markdown"].exists()
    assert summary.exists()
