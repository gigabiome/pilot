"""The uniform model API every baseline implements.

The HARNESS never imports this: it only reads prediction files (see harness.py), so
external teams can use any language or framework. This base class just makes our own
baselines consistent and removes boilerplate:

    fit(train_ids, val_ids)     learn from public labels only
    predict(ids)                -> log2 MIC + P(resistant) per isolate x antibiotic
    explain(ids, top_k)         -> ranked gene families per isolate x antibiotic (T7)

Calibration: every model predicts a log2 MIC (mu). We estimate how wrong it typically
is on the VALIDATION set (sigma, per antibiotic) and turn that into a probability:
    P(R) = P(true log2 MIC > breakpoint) = 1 - Phi((breakpoint + 0.5 - mu) / sigma)
That is what makes ECE (calibration error) meaningful for regression models.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd

from amrbench.config import Config
from amrbench.features import FeatureStore


class BenchData:
    """What a model is allowed to see: public metadata, public labels, features."""

    def __init__(self, cfg: Config):
        pub = cfg.path("release_dir") / "public"
        self.cfg = cfg
        self.meta = pd.read_csv(pub / "metadata.tsv", sep="\t").set_index("isolate_id")
        self.labels = pd.read_csv(pub / "labels_trainval.tsv", sep="\t")
        self.antibiotics: list[str] = cfg["antibiotics"]
        self.features = FeatureStore(cfg)
        self._y = self.labels.pivot(index="isolate_id", columns="antibiotic", values="log2_mic_point")

    def y(self, ids: list[str]) -> pd.DataFrame:
        """log2 MIC targets (NaN where untested). Raises if asked for non-public isolates."""
        unknown = set(ids) - set(self._y.index) - set(self.meta.index)
        if unknown:
            raise KeyError(f"unknown isolates {list(unknown)[:3]}")
        return self._y.reindex(index=ids, columns=self.antibiotics)

    def breakpoints(self, ids: list[str]) -> pd.DataFrame:
        sp = self.meta.loc[ids, "species"]
        return pd.DataFrame({ab: [self.cfg.breakpoint_log2(s, ab) for s in sp] for ab in self.antibiotics}, index=ids)


def _norm_sf(z: np.ndarray) -> np.ndarray:
    return 0.5 * np.vectorize(math.erfc)(z / math.sqrt(2))


class BaselineModel(ABC):
    name = "base"
    description = ""

    def __init__(self, data: BenchData, params: dict):
        self.data, self.params = data, params
        self.sigma: dict[str, float] = {}

    # --- to implement
    @abstractmethod
    def _fit(self, train_ids: list[str], val_ids: list[str]) -> None: ...

    @abstractmethod
    def _predict_log2(self, ids: list[str]) -> pd.DataFrame:
        """DataFrame (ids x antibiotics) of predicted log2 MIC."""

    @abstractmethod
    def _gene_scores(self, iid: str, antibiotic: str) -> dict[str, float]:
        """Gene family -> contribution to raising the MIC for this isolate."""

    # --- shared logic
    def fit(self, train_ids: list[str], val_ids: list[str]) -> "BaselineModel":
        self._fit(train_ids, val_ids)
        resid = self._predict_log2(val_ids) - self.data.y(val_ids)
        for ab in self.data.antibiotics:
            r = resid[ab].dropna()
            self.sigma[ab] = float(max(np.sqrt(np.mean(r ** 2)) if len(r) else 1.0, 0.3))
        return self

    def predict(self, ids: list[str]) -> pd.DataFrame:
        mu = self._predict_log2(ids)
        bp = self.data.breakpoints(ids)
        rows = []
        for ab in self.data.antibiotics:
            p = _norm_sf((bp[ab].to_numpy() + 0.5 - mu[ab].to_numpy()) / self.sigma[ab])
            rows.append(pd.DataFrame({"isolate_id": ids, "antibiotic": ab,
                                      "log2_mic_pred": mu[ab].to_numpy().round(3), "p_resistant": p.round(4)}))
        return pd.concat(rows, ignore_index=True)

    def explain(self, ids: list[str], top_k: int = 10) -> pd.DataFrame:
        rows = []
        for iid in ids:
            for ab in self.data.antibiotics:
                scores = self._gene_scores(iid, ab)
                ranked = sorted(((s, g) for g, s in scores.items() if s > 0), reverse=True)[:top_k]
                rows += [(iid, ab, r + 1, g, round(float(s), 5)) for r, (s, g) in enumerate(ranked)]
        return pd.DataFrame(rows, columns=["isolate_id", "antibiotic", "rank", "feature", "score"])
