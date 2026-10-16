"""Benchmark metrics. Each takes plain arrays so it can be unit-tested in isolation."""

from __future__ import annotations

import numpy as np


def essential_agreement(pred_log2, true_log2, censor) -> np.ndarray:
    """Per-pair 'within +-1 two-fold dilution' (the Appendix Q9 primary MIC metric).

    The prediction is first rounded to the nearest dilution. Censored lab values:
      '<=v' (true <= v):      agree if prediction <= v + 1
      '>v'  (true >= v + 1):  agree if prediction >= v
    """
    r = np.round(np.asarray(pred_log2, dtype=float))
    v = np.asarray(true_log2, dtype=float)
    c = np.asarray(censor)
    out = np.abs(r - v) <= 1
    out = np.where(c == "<=", r <= v + 1, out)
    out = np.where(c == ">", r >= v, out)
    return np.where(np.isnan(r), False, out)          # a missing prediction counts as wrong


def binary_metrics(y_true, y_pred) -> dict:
    """R = positive class. VME/ME are the standard AST error rates (ISO 20776-2):
    very major error = called S but truly R (dangerous); major error = called R but truly S."""
    t, p = np.asarray(y_true, bool), np.asarray(y_pred, bool)
    tp, tn = np.sum(t & p), np.sum(~t & ~p)
    fp, fn = np.sum(~t & p), np.sum(t & ~p)
    sens = tp / (tp + fn) if tp + fn else np.nan
    spec = tn / (tn + fp) if tn + fp else np.nan
    prec = tp / (tp + fp) if tp + fp else np.nan
    f1 = 2 * prec * sens / (prec + sens) if (prec + sens) else 0.0
    return {"n": int(len(t)), "balanced_accuracy": float(np.nanmean([sens, spec])), "f1": float(f1),
            "sensitivity": float(sens), "specificity": float(spec),
            "very_major_error_rate": float(fn / (tp + fn)) if tp + fn else np.nan,
            "major_error_rate": float(fp / (tn + fp)) if tn + fp else np.nan}


def expected_calibration_error(y_true, p, n_bins: int = 10) -> float:
    """ECE: average |observed R rate - predicted P(R)| over equal-width probability bins,
    weighted by bin size. 0 = perfectly calibrated."""
    y, p = np.asarray(y_true, float), np.clip(np.asarray(p, float), 0, 1)
    bins = np.minimum((p * n_bins).astype(int), n_bins - 1)
    ece = 0.0
    for b in range(n_bins):
        m = bins == b
        if m.any():
            ece += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(ece)


def topk_recall(ranked: list[str], relevant: set[str], k: int) -> float:
    """Fraction of the known determinants present in this isolate that appear in the top-k features."""
    return len(set(ranked[:k]) & relevant) / len(relevant) if relevant else np.nan


def bootstrap_ci(values: np.ndarray, groups: np.ndarray, n_boot: int = 500, seed: int = 0) -> tuple[float, float]:
    """95% CI of the mean of `values`, resampling whole isolates (pairs from one isolate are correlated)."""
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(groups, return_inverse=True)
    sums = np.bincount(inv, weights=values.astype(float))
    counts = np.bincount(inv)
    stats = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(uniq), len(uniq))
        stats.append(sums[pick].sum() / counts[pick].sum())
    lo, hi = np.percentile(stats, [2.5, 97.5])
    return float(lo), float(hi)
