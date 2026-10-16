"""Build the benchmark from validated inputs: labels -> clone clusters -> splits -> release.

The two ideas that matter most here:

1. LEAK-FREE SPLITS. Clinical isolates come in clonal outbreaks. If near-identical
   genomes sit on both sides of the train/test line, a model can "memorise the clone"
   and look far better than it will on a new patient. So we first group isolates into
   clone clusters by genomic similarity, then split whole clusters (StratifiedGroupKFold),
   while keeping species and collection site balanced across partitions.

2. HIDDEN TEST LABELS. Test labels go to release/private/ and are only ever read by the
   harness. Model code only sees release/public/. This is what keeps a public benchmark
   honest once external teams start submitting.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.model_selection import StratifiedGroupKFold

from amrbench.config import Config, slug
from amrbench.data import load_metadata, load_mic
from amrbench.features import jaccard_distance, kmer_matrix


# ----------------------------------------------------------------------------- labels

def make_labels(mic: pd.DataFrame, meta: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Add breakpoint, R/S call and a point estimate for training to every MIC.

    Censored values need care:  '<=v' is S only if v <= breakpoint, '>v' is R only if
    v >= breakpoint; otherwise the R/S call is genuinely unknown (left as NaN).
    """
    df = mic.merge(meta[["isolate_id", "species"]], on="isolate_id")
    df = df[df.antibiotic.isin(cfg["antibiotics"])].copy()
    df["breakpoint_log2"] = [cfg.breakpoint_log2(s, a) for s, a in zip(df.species, df.antibiotic)]
    v, b, c = df.log2_mic, df.breakpoint_log2, df.censor
    res = np.where(c == "=", v > b, np.nan)
    res = np.where((c == "<=") & (v <= b), 0.0, res)
    res = np.where((c == ">") & (v >= b), 1.0, res)
    df["resistant"] = res.astype(float)
    # Point estimate used as a regression target (simple choice; censored-likelihood
    # models are a possible improvement for the real data).
    df["log2_mic_point"] = np.where(c == ">", v + 1, v)
    return df[["isolate_id", "species", "antibiotic", "mic", "log2_mic", "censor",
               "log2_mic_point", "breakpoint_log2", "resistant"]]


# ----------------------------------------------------------------------------- clones

def clone_clusters(meta: pd.DataFrame, D: np.ndarray, threshold: float) -> pd.Series:
    """Average-linkage clustering of genomic distances within each species."""
    labels = pd.Series(index=meta.isolate_id, dtype=object)
    pos = {iid: i for i, iid in enumerate(meta.isolate_id)}
    for species, grp in meta.groupby("species"):
        idx = [pos[i] for i in grp.isolate_id]
        if len(idx) == 1:
            cl = np.array([0])
        else:
            cl = AgglomerativeClustering(n_clusters=None, metric="precomputed", linkage="average",
                                         distance_threshold=threshold).fit_predict(D[np.ix_(idx, idx)])
        labels.loc[grp.isolate_id] = [f"{slug(species)}_C{c:03d}" for c in cl]
    return labels


# ----------------------------------------------------------------------------- splits

def _group_fold(df: pd.DataFrame, n_splits: int, strat_cols: list[str], seed: int) -> np.ndarray:
    """Boolean mask of one fold, keeping clone clusters whole and strata balanced."""
    y = df[strat_cols].astype(str).agg("|".join, axis=1)
    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    _, test_idx = next(sgkf.split(df, y, groups=df.clone_cluster))
    mask = np.zeros(len(df), dtype=bool)
    mask[test_idx] = True
    return mask


def main_split(meta: pd.DataFrame, cfg: Config) -> dict[str, list[str]]:
    s = cfg["splits"]
    test = _group_fold(meta, round(1 / s["test_fraction"]), s["stratify_by"], cfg["seed"])
    rest = meta[~test]
    val = _group_fold(rest, round((1 - s["test_fraction"]) / s["val_fraction"]), s["stratify_by"], cfg["seed"] + 1)
    return {
        "train": sorted(rest.isolate_id[~val]),
        "val": sorted(rest.isolate_id[val]),
        "test": sorted(meta.isolate_id[test]),
    }


def random_split(meta: pd.DataFrame, cfg: Config) -> dict[str, list[str]]:
    """NAIVE split (ignores clones). Internal only - used to demonstrate leakage."""
    rng = np.random.default_rng(cfg["seed"])
    parts = {"train": [], "val": [], "test": []}
    s = cfg["splits"]
    for _, grp in meta.groupby("species"):
        ids = rng.permutation(grp.isolate_id.to_numpy())
        n_test, n_val = round(len(ids) * s["test_fraction"]), round(len(ids) * s["val_fraction"])
        parts["test"] += list(ids[:n_test])
        parts["val"] += list(ids[n_test:n_test + n_val])
        parts["train"] += list(ids[n_test + n_val:])
    return {k: sorted(v) for k, v in parts.items()}


def cross_species_split(main: dict, meta: pd.DataFrame, species: str) -> dict[str, list[str]]:
    """T6: drop one species from train/val; test = that species' isolates in the MAIN test set.

    Building it inside the main split means (a) no test labels leak into training, and
    (b) in-distribution vs out-of-distribution are compared on the very same isolates.
    """
    sp_of = meta.set_index("isolate_id").species
    return {
        "train": [i for i in main["train"] if sp_of[i] != species],
        "val": [i for i in main["val"] if sp_of[i] != species],
        "test": [i for i in main["test"] if sp_of[i] == species],
    }


# ----------------------------------------------------------------------------- audit

def audit_split(split: dict, meta: pd.DataFrame, labels: pd.DataFrame, D: np.ndarray, threshold: float) -> dict:
    """Quantify leakage and balance. This report should be published with the splits."""
    pos = {iid: i for i, iid in enumerate(meta.isolate_id)}
    m = meta.set_index("isolate_id")
    tr, te = split["train"], split["test"]
    shared_clusters = set(m.loc[tr, "clone_cluster"]) & set(m.loc[te, "clone_cluster"])
    nearest = []
    for sp_name, grp in m.loc[te].groupby("species"):
        tr_same = [i for i in tr if m.at[i, "species"] == sp_name]
        if tr_same:
            sub = D[np.ix_([pos[i] for i in grp.index], [pos[i] for i in tr_same])]
            nearest.extend(sub.min(axis=1))
    nearest = np.array(nearest)
    out = {
        "n": {k: len(v) for k, v in split.items()},
        "clone_clusters_shared_train_test": len(shared_clusters),
        "test_isolates_with_near_identical_train_genome": int((nearest < threshold).sum()) if len(nearest) else 0,
        "median_nearest_train_distance": float(np.median(nearest)) if len(nearest) else None,
        "species_counts": {k: m.loc[v, "species"].value_counts().to_dict() for k, v in split.items()},
        "site_counts": {k: m.loc[v, "collection_site"].value_counts().to_dict() for k, v in split.items()},
    }
    lab = labels.assign(part=labels.isolate_id.map({i: k for k, v in split.items() for i in v}))
    out["resistance_rate"] = (lab.dropna(subset=["part"]).groupby(["part", "antibiotic"]).resistant.mean()
                              .round(3).unstack().to_dict(orient="index"))
    return out


# ----------------------------------------------------------------------------- release

def _sha256(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_manifest(path, name: str, split: dict, audit: dict, cfg: Config, method: str, label_file) -> None:
    manifest = {
        "benchmark": cfg["benchmark_name"], "version": cfg["version"], "split": name,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "method": method, "seed": cfg["seed"], "params": cfg["splits"],
        "labels_sha256": _sha256(label_file),
        "partitions": split, "audit": audit,
    }
    path.write_text(json.dumps(manifest, indent=2, default=str))


def build(cfg: Config) -> dict:
    root, rel = cfg.path("data_root"), cfg.path("release_dir")
    if rel.exists():
        shutil.rmtree(rel)
    (rel / "public" / "splits").mkdir(parents=True)
    (rel / "private").mkdir(parents=True)

    meta = load_metadata(root)
    meta = meta[meta.species.isin(cfg["species"])]
    if cfg["qc"]["require_qc_pass"]:
        meta = meta[meta.qc_pass]
    meta = meta.reset_index(drop=True)
    labels = make_labels(load_mic(root), meta, cfg)

    s = cfg["splits"]
    sketch = kmer_matrix(meta.isolate_id.tolist(), root, s["sketch_k"], s["sketch_bits"])
    D = jaccard_distance(sketch)
    cfg.path("work_dir").mkdir(exist_ok=True)
    np.save(cfg.path("work_dir") / "distances.npy", D)
    meta["clone_cluster"] = clone_clusters(meta, D, s["clone_distance_threshold"]).loc[meta.isolate_id].to_numpy()

    main = main_split(meta, cfg)
    splits = {"main": (main, "clone-grouped, stratified by " + "+".join(s["stratify_by"]))}
    for species in s["cross_species_holdout"]:
        splits[f"T6_{slug(species)}"] = (cross_species_split(main, meta, species), f"main split without {species}")

    # Public: metadata + train/val labels. Private: everything (the harness needs test labels).
    trainval = set(main["train"]) | set(main["val"])
    part = {i: k for k, v in main.items() for i in v}
    meta.to_csv(rel / "public" / "metadata.tsv", sep="\t", index=False)
    public_labels = labels[labels.isolate_id.isin(trainval)].assign(partition=lambda d: d.isolate_id.map(part))
    public_labels.to_csv(rel / "public" / "labels_trainval.tsv", sep="\t", index=False)
    labels.assign(partition=labels.isolate_id.map(part)).to_csv(rel / "private" / "labels_all.tsv", sep="\t", index=False)

    audits = {}
    for name, (split, method) in splits.items():
        audits[name] = audit_split(split, meta, labels, D, s["clone_distance_threshold"])
        write_manifest(rel / "public" / "splits" / f"{name}.json", name, split, audits[name], cfg, method,
                       rel / "private" / "labels_all.tsv")
    # Internal leakage demo split (not part of the public release).
    rnd = random_split(meta, cfg)
    audits["random_INTERNAL"] = audit_split(rnd, meta, labels, D, s["clone_distance_threshold"])
    write_manifest(rel / "private" / "random_split_INTERNAL.json", "random_INTERNAL", rnd,
                   audits["random_INTERNAL"], cfg, "naive random (leaky, demo only)", rel / "private" / "labels_all.tsv")

    (rel / "public" / "tasks.json").write_text(json.dumps(
        {"tasks": cfg["tasks"], "antibiotics": cfg["antibiotics"], "breakpoints_mg_per_L": cfg["breakpoints"],
         "mic_panel": cfg["mic_panel"]}, indent=2))
    return {"n_isolates": len(meta), "n_labels": len(labels), "n_clusters": meta.clone_cluster.nunique(),
            "audits": audits}
