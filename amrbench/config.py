"""Load config/benchmark.yaml and provide small helpers used everywhere."""

from __future__ import annotations

import math
import re
from pathlib import Path

import yaml

DEFAULT_CONFIG = Path(__file__).resolve().parent.parent / "config" / "benchmark.yaml"


class Config(dict):
    """A dict with a few convenience accessors. Paths are resolved relative to the project root."""

    def __init__(self, data: dict, root: Path):
        super().__init__(data)
        self.root = root

    def path(self, key: str) -> Path:
        p = Path(self[key])
        return p if p.is_absolute() else self.root / p

    def breakpoint_log2(self, species: str, antibiotic: str) -> float:
        """log2 of the R> breakpoint (mg/L). Species-specific value wins over the default."""
        bps = self["breakpoints"]
        value = bps.get(species, {}).get(antibiotic, bps["default"][antibiotic])
        return math.log2(value)


def load_config(path: str | Path | None = None) -> Config:
    path = Path(path) if path else DEFAULT_CONFIG
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return Config(data, root=path.resolve().parent.parent)


def slug(text: str) -> str:
    """'Klebsiella pneumoniae' -> 'Klebsiella_pneumoniae' (safe for filenames)."""
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
