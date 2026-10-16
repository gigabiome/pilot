"""Command line entry point:  python -m amrbench <step>   (or `amrbench <step>` after pip install)

Steps, in order:
  generate-mock   write a synthetic dataset to data_root          (skip once real data exists)
  validate        check data_root against the data contract
  build           labels, clone clusters, leak-free splits, public/private release
  features        k-mer matrix, pan-genome matrix, gene graphs
  train           train baselines and write their submissions
  evaluate        score every submission with the harness
  leaderboard     collect scores into results/leaderboard.md
  leakage-demo    internal: how much would a naive random split have inflated scores?
  run-all         all of the above
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

from amrbench.config import load_config


def _step(title):
    print(f"\n=== {title} " + "=" * max(0, 60 - len(title)))


def cmd_generate(cfg, args):
    from amrbench.mock.generate import generate
    _step("generate-mock")
    out = generate(cfg)
    print(f"mock dataset written to {out}")


def cmd_validate(cfg, args):
    from amrbench.data import validate
    _step("validate")
    for note in validate(cfg.path("data_root"), cfg):
        print(" ", note)


def cmd_build(cfg, args):
    from amrbench.build import build
    _step("build")
    r = build(cfg)
    print(f"  {r['n_isolates']} QC-passed isolates, {r['n_labels']} MIC labels, {r['n_clusters']} clone clusters")
    for name, a in r["audits"].items():
        med = a["median_nearest_train_distance"]
        print(f"  split {name:<22} sizes={a['n']}  shared clones={a['clone_clusters_shared_train_test']}  "
              f"near-identical test genomes={a['test_isolates_with_near_identical_train_genome']}"
              + (f"  median nearest-train distance={med:.3f}" if med is not None else ""))


def cmd_features(cfg, args):
    from amrbench.features import build_features
    _step("features")
    ids = pd.read_csv(cfg.path("release_dir") / "public" / "metadata.tsv", sep="\t").isolate_id.tolist()
    t0 = time.time()
    out = build_features(cfg, ids)
    print(f"  features for {len(ids)} isolates written to {out} ({time.time() - t0:.1f}s)")


def cmd_train(cfg, args):
    from amrbench.baselines import MODELS, run_baseline
    from amrbench.models.base import BenchData
    _step("train baselines")
    data = BenchData(cfg)
    for key in (args.models or list(MODELS)):
        run_baseline(cfg, key, data)


def cmd_evaluate(cfg, args):
    from amrbench.harness import evaluate
    _step("evaluate")
    subs = [Path(args.submission)] if getattr(args, "submission", None) else \
        sorted(p for p in cfg.path("submissions_dir").iterdir() if (p / "predictions_main.tsv").exists())
    for sub in subs:
        s = evaluate(cfg, sub)
        t = s["tasks"]
        print(f"  {sub.name:<16} T1 EA={t['T1']['essential_agreement']:.3f}  T5 F1={t['T5']['f1']:.3f}  "
              f"T7 recall={t['T7'].get('top' + str(cfg['tasks']['T7'].get('top_k', 5)) + '_recall')}  -> results/{sub.name}/scores.json")


def cmd_leaderboard(cfg, args):
    from amrbench.leaderboard import build_leaderboard
    _step("leaderboard")
    df = build_leaderboard(cfg)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(df.to_string(index=False))
    print(f"\n  written to {cfg.path('results_dir') / 'leaderboard.md'}")


def cmd_leakage(cfg, args):
    from amrbench.experiments import leakage_demo
    _step("leakage-demo (internal)")
    print(leakage_demo(cfg).to_string(index=False))


def cmd_all(cfg, args):
    for fn in (cmd_generate, cmd_validate, cmd_build, cmd_features, cmd_train, cmd_evaluate, cmd_leaderboard,
               cmd_leakage):
        fn(cfg, args)


COMMANDS = {"generate-mock": cmd_generate, "validate": cmd_validate, "build": cmd_build, "features": cmd_features,
            "train": cmd_train, "evaluate": cmd_evaluate, "leaderboard": cmd_leaderboard,
            "leakage-demo": cmd_leakage, "run-all": cmd_all}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="amrbench", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=COMMANDS)
    ap.add_argument("--config", default=None, help="path to benchmark.yaml (default: config/benchmark.yaml)")
    ap.add_argument("--models", nargs="*", help="train: subset of kmer pangenome gnn")
    ap.add_argument("--submission", help="evaluate: score a single submission folder")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    COMMANDS[args.command](cfg, args)


if __name__ == "__main__":
    sys.exit(main())
