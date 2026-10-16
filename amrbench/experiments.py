"""Internal experiments (not part of the public benchmark).

leakage_demo: train the same model on a naive random split and on the clone-grouped
split, and compare test scores. The gap is the over-optimism a leaky split would have
reported. It uses PRIVATE labels, which is fine for an internal methods check but is
exactly what an external submission must never do.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from amrbench.config import Config
from amrbench.metrics import essential_agreement
from amrbench.models.base import BenchData
from amrbench.baselines import MODELS, load_split


class _PrivateData(BenchData):
    """BenchData that can see all labels - internal experiments only."""

    def __init__(self, cfg):
        super().__init__(cfg)
        self.labels = pd.read_csv(cfg.path("release_dir") / "private" / "labels_all.tsv", sep="\t")
        self._y = self.labels.pivot(index="isolate_id", columns="antibiotic", values="log2_mic_point")


def leakage_demo(cfg: Config, model_key: str = "pangenome") -> pd.DataFrame:
    data = _PrivateData(cfg)
    rel = cfg.path("release_dir")
    splits = {
        "random (leaky)": json.loads((rel / "private" / "random_split_INTERNAL.json").read_text())["partitions"],
        "clone-grouped": load_split(cfg, "main"),
    }
    rows = []
    for name, s in splits.items():
        m = MODELS[model_key](data, cfg["models"][model_key]).fit(s["train"], s["val"])
        pred = m.predict(s["test"])
        df = data.labels.merge(pred, on=["isolate_id", "antibiotic"])
        ea = essential_agreement(df.log2_mic_pred, df.log2_mic, df.censor)
        rows.append({"split": name, "model": m.name, "test_isolates": len(s["test"]),
                     "essential_agreement": round(float(np.mean(ea)), 4)})
    out = pd.DataFrame(rows)
    out.to_csv(cfg.path("results_dir") / "leakage_demo.tsv", sep="\t", index=False)
    return out
