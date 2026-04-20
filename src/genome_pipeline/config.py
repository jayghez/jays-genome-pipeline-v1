"""YAML configuration loading and objective resolution."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from .io_utils import read_gene_list


DEFAULT_PIPELINE_CONFIG = Path("configs/pipeline.yaml")
DEFAULT_OBJECTIVES_CONFIG = Path("configs/objectives.yaml")


def load_yaml(path: str | Path) -> dict[str, Any]:
    config_path = Path(path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    data = yaml.safe_load(config_path.read_text()) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Expected YAML mapping in {config_path}")
    return data


def load_pipeline_config(path: str | Path = DEFAULT_PIPELINE_CONFIG) -> dict[str, Any]:
    return load_yaml(path)


def load_objectives(path: str | Path = DEFAULT_OBJECTIVES_CONFIG) -> dict[str, dict[str, Any]]:
    objectives_path = Path(path)
    data = load_yaml(objectives_path)
    raw_objectives = data.get("objectives", {})
    if not isinstance(raw_objectives, dict):
        raise ValueError("objectives.yaml must contain an 'objectives' mapping")

    resolved: dict[str, dict[str, Any]] = {}
    for key, objective in raw_objectives.items():
        if not isinstance(objective, dict):
            raise ValueError(f"Objective {key} must be a mapping")
        copy = deepcopy(objective)
        genes = [gene.upper() for gene in copy.get("genes", [])]
        gene_list_file = copy.get("gene_list_file")
        if gene_list_file:
            genes.extend(read_gene_list(objectives_path.parent / gene_list_file))
        copy["genes"] = sorted(set(genes))
        copy.setdefault("name", key)
        copy.setdefault("aliases", [])
        copy.setdefault("workflow", "disease_risk")
        resolved[key] = copy
    return resolved


def get_objective(objectives: dict[str, dict[str, Any]], name: str) -> dict[str, Any]:
    wanted = name.lower()
    for key, objective in objectives.items():
        aliases = [str(alias).lower() for alias in objective.get("aliases", [])]
        candidates = {key.lower(), str(objective.get("name", "")).lower(), *aliases}
        if wanted in candidates:
            return deepcopy(objective)
    available = ", ".join(sorted(objectives))
    raise KeyError(f"Unknown objective '{name}'. Available objectives: {available}")


def resource_path(config: dict[str, Any], key: str) -> Path | None:
    value = config.get("resources", {}).get(key)
    if not value:
        return None
    return Path(value)
