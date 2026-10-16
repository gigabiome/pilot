"""Collect results/*/scores.json into one leaderboard table (markdown + TSV)."""

from __future__ import annotations

import json

import pandas as pd

from amrbench.config import Config


def _pct(x):
    return None if x is None else round(100 * x, 1)


def summarise(s: dict, top_k: int) -> dict:
    t = s["tasks"]
    t6 = [v for v in t["T6"].values() if isinstance(v, dict) and "ea_drop_pp" in v]
    return {
        "model": s["submission"].get("model", "?"),
        "T1 EA %": _pct(t["T1"]["essential_agreement"]),
        "T2 EA %": _pct(t["T2"]["essential_agreement"]),
        "T3 EA %": _pct(t["T3"]["essential_agreement"]),
        "T4 min-species EA %": _pct(t["T4"].get("min_species_ea")),
        "T5 F1": t["T5"]["f1"],
        "T5 sens": t["T5"]["sensitivity"],
        "T5 VME %": _pct(t["T5"]["very_major_error_rate"]),
        "T6 EA drop pp": round(sum(v["ea_drop_pp"] for v in t6) / len(t6), 1) if t6 else None,
        "T6 ECE": round(sum(v["ece_held_out"] for v in t6) / len(t6), 4) if t6 else None,
        f"T7 top{top_k} recall": t["T7"].get(f"top{top_k}_recall"),
        "T7 recall (mutations)": t["T7"].get(f"top{top_k}_recall_by_mechanism", {}).get("point_mutation"),
        "T7 ECE": t["T7"].get("ece"),
        "decoy rate (mock)": next((v for kk, v in t["T7"].items() if kk.endswith("MOCK_ONLY")), None),
    }


def build_leaderboard(cfg: Config) -> pd.DataFrame:
    tasks = cfg["tasks"]
    k = tasks["T7"].get("top_k", 5)
    rows = [summarise(json.loads(p.read_text()), k) for p in sorted(cfg.path("results_dir").glob("*/scores.json"))]
    target = {"model": "TARGET (App. Q9)",
              "T1 EA %": 100 * tasks["T1"]["target"]["essential_agreement"],
              "T2 EA %": 100 * tasks["T2"]["target"]["essential_agreement"],
              "T3 EA %": 100 * tasks["T3"]["target"]["essential_agreement"],
              "T4 min-species EA %": 100 * tasks["T4"]["target"]["essential_agreement"],
              "T5 F1": tasks["T5"]["target"]["f1"], "T5 sens": tasks["T5"]["target"]["sensitivity"],
              "T6 EA drop pp": f"<= {tasks['T6']['target']['ea_drop_pp_max']}",
              "T6 ECE": f"<= {tasks['T6']['target']['ece_max']}",
              f"T7 top{k} recall": tasks["T7"]["target"]["topk_recall"],
              "T7 ECE": f"<= {tasks['T7']['target']['ece_max']}"}
    df = pd.DataFrame(rows + [target])
    out = cfg.path("results_dir")
    df.to_csv(out / "leaderboard.tsv", sep="\t", index=False)
    md = ["# AMRBench pilot leaderboard (MOCK DATA - numbers are not scientific results)", "",
          "| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    for _, r in df.iterrows():
        md.append("| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |")
    (out / "leaderboard.md").write_text("\n".join(md) + "\n")
    return df
