"""The evaluation harness: score a SUBMISSION (files on disk) against the hidden test labels.

Submission format (the public contract for external teams)
-----------------------------------------------------------
submissions/<name>/
  submission.json              {"model": ..., "team": ..., "description": ...}
  predictions_main.tsv         isolate_id, antibiotic, log2_mic_pred, p_resistant
                               for every test isolate x antibiotic in splits/main.json
  attributions_main.tsv        (optional, T7) isolate_id, antibiotic, rank, feature
  predictions_T6_<species>.tsv (optional, T6) same columns, model trained on splits/T6_<species>.json

`p_resistant` is optional; without it T5 uses log2_mic_pred vs breakpoint and ECE is skipped.
Missing predictions count as wrong (coverage is reported).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from amrbench.config import Config, slug
from amrbench.data import load_amrfinder, load_known_determinants
from amrbench.metrics import (binary_metrics, bootstrap_ci, essential_agreement,
                              expected_calibration_error, topk_recall)

PRED_COLUMNS = ["isolate_id", "antibiotic", "log2_mic_pred"]


class SubmissionError(ValueError):
    pass


def load_predictions(path: Path, test_ids: set[str], antibiotics: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "log2_mic_pred" not in df and "mic_pred" in df:
        df["log2_mic_pred"] = np.log2(df["mic_pred"].astype(float))   # mg/L also accepted
    missing = [c for c in PRED_COLUMNS if c not in df.columns]
    if missing:
        raise SubmissionError(f"{path.name}: missing columns {missing}")
    if df.duplicated(["isolate_id", "antibiotic"]).any():
        raise SubmissionError(f"{path.name}: duplicate isolate_id x antibiotic rows")
    extra = set(df.isolate_id) - test_ids
    if extra:
        raise SubmissionError(f"{path.name}: {len(extra)} isolates are not in this test split")
    if "p_resistant" in df and not df.p_resistant.dropna().between(0, 1).all():
        raise SubmissionError(f"{path.name}: p_resistant must be in [0, 1]")
    return df[df.antibiotic.isin(antibiotics)]


def _join(labels: pd.DataFrame, pred: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in ["isolate_id", "antibiotic", "log2_mic_pred", "p_resistant"] if c in pred]
    df = labels.merge(pred[cols], on=["isolate_id", "antibiotic"], how="left")
    df["ea"] = essential_agreement(df.log2_mic_pred, df.log2_mic, df.censor)
    return df


def _ea_summary(df: pd.DataFrame) -> dict:
    if df.empty:
        return {"essential_agreement": None, "n": 0}
    lo, hi = bootstrap_ci(df.ea.to_numpy(), df.isolate_id.to_numpy())
    return {"essential_agreement": round(float(df.ea.mean()), 4), "ci95": [round(lo, 4), round(hi, 4)],
            "n_pairs": int(len(df)), "n_isolates": int(df.isolate_id.nunique()),
            "coverage": round(float(df.log2_mic_pred.notna().mean()), 4),
            "per_antibiotic": df.groupby("antibiotic").ea.mean().round(4).to_dict()}


def _resistance_call(df: pd.DataFrame) -> np.ndarray:
    if "p_resistant" in df and df.p_resistant.notna().all():
        return (df.p_resistant >= 0.5).to_numpy()
    return (df.log2_mic_pred.fillna(-np.inf) > df.breakpoint_log2).to_numpy()


def _ece(df: pd.DataFrame) -> float | None:
    d = df.dropna(subset=["resistant"])
    if "p_resistant" not in d or d.p_resistant.isna().any() or d.empty:
        return None
    return round(expected_calibration_error(d.resistant, d.p_resistant), 4)


# ----------------------------------------------------------------------------- task scorers

def score_mic(task, df):
    d = df[df.species.isin(task["species"])]
    out = _ea_summary(d)
    if task.get("per_species"):
        out["per_species"] = {s: _ea_summary(g) for s, g in d.groupby("species")}
        out["min_species_ea"] = min(v["essential_agreement"] for v in out["per_species"].values())
    return out


def score_binary(task, df):
    d = df.dropna(subset=["resistant"])
    out = binary_metrics(d.resistant.astype(bool), _resistance_call(d))
    out["ece"] = _ece(d)
    out["per_antibiotic_f1"] = {ab: round(binary_metrics(g.resistant.astype(bool), _resistance_call(g))["f1"], 4)
                                for ab, g in d.groupby("antibiotic")}
    return {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out.items()}


def score_cross_species(task, df_main, sub_dir, cfg, labels, splits_dir):
    out = {}
    for species in cfg["splits"]["cross_species_holdout"]:
        name = f"T6_{slug(species)}"
        path = sub_dir / f"predictions_{name}.tsv"
        test = set(json.loads((splits_dir / f"{name}.json").read_text())["partitions"]["test"])
        if not path.exists():
            out[species] = {"status": "not submitted"}
            continue
        ood = _join(labels[labels.isolate_id.isin(test)], load_predictions(path, test, cfg["antibiotics"]))
        ind = df_main[df_main.isolate_id.isin(test)]
        ea_in, ea_out = ind.ea.mean(), ood.ea.mean()
        out[species] = {"ea_in_distribution": round(float(ea_in), 4), "ea_held_out": round(float(ea_out), 4),
                        "ea_drop_pp": round(100 * float(ea_in - ea_out), 2), "ece_held_out": _ece(ood),
                        "n_pairs": int(len(ood))}
    return out


def score_attribution(task, df, sub_dir, cfg):
    path = sub_dir / "attributions_main.tsv"
    k = task.get("top_k", 5)
    out = {"ece": _ece(df)}
    if not path.exists():
        return {**out, "status": "no attributions submitted"}
    attr = pd.read_csv(path, sep="\t").sort_values("rank")
    ranked = attr.groupby(["isolate_id", "antibiotic"]).feature.apply(list).to_dict()

    root = cfg.path("data_root")
    known = load_known_determinants(root)
    det_by_ab = known.groupby("antibiotic").gene_family.apply(set).to_dict()
    mechanism = dict(zip(known.gene_family, known.mechanism))
    carried = load_amrfinder(root).groupby("isolate_id").gene_family.apply(set).to_dict()

    # Only score resistant predictions that a known determinant can explain
    # (intrinsic resistance has no gene to point at).
    recalls, decoy_hits, by_mech = [], [], {}
    truth_file = root / "truth" / "ground_truth.json"
    decoys = set(json.loads(truth_file.read_text())["decoys_plasmid_backbone"]) if truth_file.exists() else None
    for r in df[df.resistant == 1].itertuples():
        relevant = det_by_ab.get(r.antibiotic, set()) & carried.get(r.isolate_id, set())
        if not relevant:
            continue
        top = ranked.get((r.isolate_id, r.antibiotic), [])
        recalls.append(topk_recall(top, relevant, k))
        for det in relevant:                      # found / not found, per mechanism type
            by_mech.setdefault(mechanism[det], []).append(det in top[:k])
        if decoys is not None:
            decoy_hits.append(len(set(top[:k]) & decoys) / k)
    out.update({f"top{k}_recall": round(float(np.mean(recalls)), 4) if recalls else None, "n_pairs": len(recalls),
                f"top{k}_recall_by_mechanism": {m: round(float(np.mean(v)), 4) for m, v in by_mech.items()}})
    if decoys is not None:
        out[f"top{k}_decoy_rate_MOCK_ONLY"] = round(float(np.mean(decoy_hits)), 4) if decoy_hits else None
    return out


# ----------------------------------------------------------------------------- entry point

def evaluate(cfg: Config, sub_dir: Path) -> dict:
    rel = cfg.path("release_dir")
    splits_dir = rel / "public" / "splits"
    labels = pd.read_csv(rel / "private" / "labels_all.tsv", sep="\t")
    test = set(json.loads((splits_dir / "main.json").read_text())["partitions"]["test"])
    test_labels = labels[labels.isolate_id.isin(test)]
    df = _join(test_labels, load_predictions(sub_dir / "predictions_main.tsv", test, cfg["antibiotics"]))

    meta = json.loads((sub_dir / "submission.json").read_text()) if (sub_dir / "submission.json").exists() else {}
    scores = {"submission": meta, "benchmark": cfg["benchmark_name"], "version": cfg["version"], "tasks": {}}
    for tid, task in cfg["tasks"].items():
        kind = task["kind"]
        if kind == "mic":
            s = score_mic(task, df)
        elif kind == "binary":
            s = score_binary(task, df)
        elif kind == "cross_species":
            s = score_cross_species(task, df, sub_dir, cfg, labels, splits_dir)
        elif kind == "attribution":
            s = score_attribution(task, df, sub_dir, cfg)
        else:
            raise ValueError(f"unknown task kind {kind}")
        scores["tasks"][tid] = {"name": task["name"], "target": task.get("target"), **s}

    out = cfg.path("results_dir") / sub_dir.name
    out.mkdir(parents=True, exist_ok=True)
    (out / "scores.json").write_text(json.dumps(scores, indent=2, default=float))
    df.to_csv(out / "test_pairs_scored.tsv", sep="\t", index=False)
    return scores
