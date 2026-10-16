"""Train the three named baselines and write their submissions.

Our baselines go through exactly the same door as an external team would: they write
prediction files, and the harness scores those files. Nothing is special-cased.
"""

from __future__ import annotations

import json
import time

from amrbench.config import Config, slug
from amrbench.models.base import BenchData
from amrbench.models.gnn import GNNModel
from amrbench.models.kmer import KmerModel
from amrbench.models.pangenome import PangenomeModel

MODELS = {"kmer": KmerModel, "pangenome": PangenomeModel, "gnn": GNNModel}


def load_split(cfg: Config, name: str) -> dict:
    return json.loads((cfg.path("release_dir") / "public" / "splits" / f"{name}.json").read_text())["partitions"]


def run_baseline(cfg: Config, key: str, data: BenchData | None = None, log=print) -> None:
    data = data or BenchData(cfg)
    cls = MODELS[key]
    out = cfg.path("submissions_dir") / cls.name
    out.mkdir(parents=True, exist_ok=True)
    top_k = max(t.get("top_k", 5) for t in cfg["tasks"].values()) * 2

    main = load_split(cfg, "main")
    t0 = time.time()
    model = cls(data, cfg["models"][key]).fit(main["train"], main["val"])
    model.predict(main["test"]).to_csv(out / "predictions_main.tsv", sep="\t", index=False)
    model.explain(main["test"], top_k=top_k).to_csv(out / "attributions_main.tsv", sep="\t", index=False)
    log(f"  {cls.name}: main split trained + predicted in {time.time() - t0:.1f}s")

    for species in cfg["splits"]["cross_species_holdout"]:
        name = f"T6_{slug(species)}"
        s = load_split(cfg, name)
        t0 = time.time()
        m = cls(data, cfg["models"][key]).fit(s["train"], s["val"])
        m.predict(s["test"]).to_csv(out / f"predictions_{name}.tsv", sep="\t", index=False)
        log(f"  {cls.name}: {name} trained + predicted in {time.time() - t0:.1f}s")

    (out / "submission.json").write_text(json.dumps({
        "model": cls.name, "team": "AMRBench baselines", "description": cls.description,
        "params": cfg["models"][key]}, indent=2))
