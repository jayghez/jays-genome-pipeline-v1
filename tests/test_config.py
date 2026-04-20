from genome_pipeline.config import get_objective, load_objectives


def test_load_objectives_resolves_gene_files():
    objectives = load_objectives("configs/objectives.yaml")
    secondary = get_objective(objectives, "secondary_findings")
    pgx = get_objective(objectives, "pharmacogenomics")

    assert "BRCA1" in secondary["genes"]
    assert "CYP2D6" in pgx["genes"]
    assert pgx["workflow"] == "pharmacogenomics"
