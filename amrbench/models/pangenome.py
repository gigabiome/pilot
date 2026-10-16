"""Baseline 2 - pan-genome model.

Input : gene family presence/absence (Panaroo matrix), plus species one-hot.
Model : one random forest per antibiotic on log2 MIC.
Sees  : gene gain/loss and interactions between genes. BLIND to point mutations,
        because a mutated gyrA is still "gyrA present" - an expected, instructive
        weakness (real fix: add Snippy/AMRFinder variant calls as extra columns).
Explain: per isolate, impurity importance of each gene it carries (direction-agnostic;
        a SHAP explainer is the natural upgrade for the real benchmark).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

from amrbench.models.base import BaselineModel


class PangenomeModel(BaselineModel):
    name = "pangenome_rf"
    description = "Random forest on gene presence/absence (Panaroo-style matrix)"

    def _matrix(self, ids):
        X = self.data.features.pangenome(ids)[self.genes]
        sp = pd.get_dummies(self.data.meta.loc[ids, "species"]).reindex(columns=self.species, fill_value=0)
        return np.hstack([X.to_numpy(np.float32), sp.to_numpy(np.float32)])

    def _fit(self, train_ids, val_ids):
        X = self.data.features.pangenome(train_ids)
        prev = X.mean(axis=0)
        lo = self.params["min_prevalence"]
        self.genes = list(prev.index[(prev >= lo) & (prev <= 1 - lo)])
        self.species = sorted(self.data.meta.species.unique())
        Xtr, ytr = self._matrix(train_ids), self.data.y(train_ids)
        self.models = {}
        for ab in self.data.antibiotics:
            m = ytr[ab].notna().to_numpy()
            self.models[ab] = RandomForestRegressor(
                n_estimators=self.params["n_estimators"], min_samples_leaf=self.params["min_samples_leaf"],
                max_features=self.params["max_features"], n_jobs=-1, random_state=0,
            ).fit(Xtr[m], ytr[ab][m])

    def _predict_log2(self, ids):
        X = self._matrix(ids)
        return pd.DataFrame({ab: m.predict(X) for ab, m in self.models.items()}, index=ids)

    def _gene_scores(self, iid, antibiotic):
        imp = self.models[antibiotic].feature_importances_[:len(self.genes)]
        present = self.data.features.pangenome([iid])[self.genes].to_numpy()[0] > 0
        return {g: float(s) for g, s, p in zip(self.genes, imp, present) if p}
