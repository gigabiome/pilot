"""Unit tests for the pieces whose correctness the whole benchmark depends on."""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from sklearn.linear_model import Ridge

from amrbench.build import make_labels
from amrbench.config import load_config
from amrbench.data import parse_mic
from amrbench.features import jaccard_distance, kmer_bins
from amrbench.metrics import binary_metrics, essential_agreement, expected_calibration_error, topk_recall
from amrbench.models.kmer import _ridge_dual


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


# --------------------------------------------------------------------------- MIC parsing

@pytest.mark.parametrize("raw,expected", [
    ("0.25", (-2.0, "=")), ("<=0.0625", (-4.0, "<=")), (">64", (6.0, ">")),
    ("<0.125", (-4.0, "<=")), (">=128", (6.0, ">")), ("≤0.5", (-1.0, "<=")),
])
def test_parse_mic(raw, expected):
    assert parse_mic(raw) == expected


def test_parse_mic_rejects_non_dilution():
    with pytest.raises(ValueError):
        parse_mic("3")


# --------------------------------------------------------------------------- labels

def test_censored_resistance_calls():
    cfg = load_config()
    bp = 2 ** cfg.breakpoint_log2("Escherichia coli", "meropenem")          # 8 mg/L -> log2 3
    meta = pd.DataFrame({"isolate_id": list("abcde"), "species": "Escherichia coli"})
    mic = pd.DataFrame({"isolate_id": list("abcde"), "antibiotic": "meropenem", "mic": ["8", "16", "<=8", "<=16", ">4"],
                        "log2_mic": [3.0, 4.0, 3.0, 4.0, 2.0], "censor": ["=", "=", "<=", "<=", ">"]})
    lab = make_labels(mic, meta, cfg).set_index("isolate_id")
    assert bp == 8
    assert lab.resistant["a"] == 0 and lab.resistant["b"] == 1           # R if MIC > breakpoint
    assert lab.resistant["c"] == 0                                         # <=8 is surely S
    assert np.isnan(lab.resistant["d"]) and np.isnan(lab.resistant["e"])   # genuinely ambiguous


# --------------------------------------------------------------------------- metrics

def test_essential_agreement_with_censoring():
    pred = [3.4, 5.0, -10, 2.0, 6.0, np.nan]
    true = [3, 3, -4, -4, 6, 0]
    cens = ["=", "=", "<=", "<=", ">", "="]
    assert essential_agreement(pred, true, cens).tolist() == [True, False, True, False, True, False]


def test_binary_metrics():
    m = binary_metrics([1, 1, 0, 0], [1, 0, 0, 1])
    assert m["sensitivity"] == 0.5 and m["specificity"] == 0.5 and m["very_major_error_rate"] == 0.5


def test_ece_perfect_and_bad():
    y = np.array([0] * 50 + [1] * 50)
    assert expected_calibration_error(y, y.astype(float)) == 0
    assert expected_calibration_error(y, 1 - y.astype(float)) == pytest.approx(1.0)


def test_topk_recall():
    assert topk_recall(["a", "b", "c"], {"a", "z"}, k=2) == 0.5


# --------------------------------------------------------------------------- features

def test_kmers_are_strand_independent():
    rng = np.random.default_rng(0)
    seq = "".join(rng.choice(list("ACGT"), 500))
    assert np.array_equal(kmer_bins(seq, 15, 18), kmer_bins(revcomp(seq), 15, 18))


def test_single_snp_changes_k_kmers():
    seq = "".join(np.random.default_rng(1).choice(list("ACGT"), 400))
    mut = seq[:200] + ("A" if seq[200] != "A" else "C") + seq[201:]
    a, b = set(kmer_bins(seq, 15, 30)), set(kmer_bins(mut, 15, 30))
    assert len(b - a) == 15


def test_jaccard_distance():
    X = sp.csr_matrix(np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 1, 0]], dtype=np.float32))
    D = jaccard_distance(X)
    assert D[0, 1] == 0 and D[0, 2] == 1 and D[0, 3] == pytest.approx(2 / 3)


# --------------------------------------------------------------------------- models

def test_dual_ridge_matches_sklearn():
    rng = np.random.default_rng(0)
    X = sp.csr_matrix((rng.random((40, 300)) < 0.3).astype(np.float64))
    y = rng.normal(size=40)
    w, b = _ridge_dual(X, y, 5.0)
    ref = Ridge(alpha=5.0).fit(X.toarray(), y)
    assert np.allclose(w, ref.coef_, atol=1e-6) and b == pytest.approx(ref.intercept_, abs=1e-6)
