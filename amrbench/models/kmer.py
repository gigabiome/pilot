"""Baseline 1 - k-mer model.

Input : presence/absence of every hashed canonical k-mer in the assembly.
Model : one ridge regression per antibiotic on log2 MIC (alpha tuned on validation),
        solved in its dual form because genomes << k-mers.
Sees  : acquired genes AND point mutations (a SNP creates k new k-mers).
Explain: per isolate, sum the coefficients of the k-mers inside each of its genes.

Real-data upgrades to consider: unitigs instead of hashed k-mers, sparse models
(elastic net) or gradient boosting, and population-structure correction.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from amrbench.models.base import BaselineModel


def _ridge_dual(X, y, alpha, K=None):
    """Ridge regression with intercept, solved in the n x n 'kernel' space.

    Same answer as sklearn's Ridge, but with n_genomes << n_kmers it is much faster:
    w = Xc^T (Xc Xc^T + alpha I)^-1 (y - mean(y)),  with Xc = column-centred X.
    K = X X^T can be passed in so it is computed once and reused.
    """
    mu = np.asarray(X.mean(axis=0)).ravel()
    K = (X @ X.T).toarray() if K is None else K
    H = np.eye(len(y)) - 1.0 / len(y)
    Kc = H @ K @ H                                             # = Xc Xc^T without densifying X
    a = np.linalg.solve(Kc + alpha * np.eye(len(y)), y - y.mean())
    w = X.T @ a - mu * a.sum()
    return w, float(y.mean() - mu @ w)


class KmerModel(BaselineModel):
    name = "kmer_ridge"
    description = "Ridge regression on hashed canonical k-mer presence"

    def _fit(self, train_ids, val_ids):
        X = self.data.features.kmer(train_ids)
        prev = np.asarray(X.mean(axis=0)).ravel()
        lo = self.params["min_prevalence"]
        self.keep = np.flatnonzero((prev >= lo) & (prev <= 1 - lo))      # drop rare/ubiquitous k-mers
        Xtr, Xva = X[:, self.keep], self.data.features.kmer(val_ids)[:, self.keep]
        ytr, yva = self.data.y(train_ids), self.data.y(val_ids)
        K = (Xtr @ Xtr.T).toarray()                               # genome x genome similarity
        self.coef, self.intercept = {}, {}
        for ab in self.data.antibiotics:
            mtr, mva = ytr[ab].notna().to_numpy(), yva[ab].notna().to_numpy()
            Kab = K[np.ix_(mtr, mtr)]
            best = None
            for alpha in self.params["alphas"]:
                w, b = _ridge_dual(Xtr[mtr], ytr[ab][mtr].to_numpy(), alpha, Kab)
                err = np.mean(np.abs(Xva[mva] @ w + b - yva[ab][mva].to_numpy()))
                if best is None or err < best[0]:
                    best = (err, alpha)
            w, b = _ridge_dual(Xtr[mtr], ytr[ab][mtr].to_numpy(), best[1], Kab)
            self.coef[ab] = np.zeros(X.shape[1])
            self.coef[ab][self.keep] = w
            self.intercept[ab] = b

    def _predict_log2(self, ids):
        X = self.data.features.kmer(ids)
        return pd.DataFrame({ab: X @ self.coef[ab] + self.intercept[ab] for ab in self.coef}, index=ids)

    def _gene_scores(self, iid, antibiotic):
        coef = self.coef[antibiotic]
        scores: dict[str, float] = {}
        for fam, bins in self.data.features.gene_kmer_bins(iid):
            scores[fam] = max(scores.get(fam, -np.inf), float(coef[bins].sum()))
        return scores
