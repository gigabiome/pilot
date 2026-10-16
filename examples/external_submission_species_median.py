"""Example of what an EXTERNAL team does - using only the public release files.

This "model" ignores the genome entirely: it predicts, for each species x antibiotic,
the median MIC seen in training. It is the floor every real model must beat
(if a genome model can't beat it, the genome isn't being used).

    python examples/external_submission_species_median.py
    python -m amrbench evaluate --submission submissions/species_median
"""

import json
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PUBLIC = ROOT / "release" / "public"
OUT = ROOT / "submissions" / "species_median"


def main():
    meta = pd.read_csv(PUBLIC / "metadata.tsv", sep="\t").set_index("isolate_id")
    labels = pd.read_csv(PUBLIC / "labels_trainval.tsv", sep="\t")
    tasks = json.loads((PUBLIC / "tasks.json").read_text())
    OUT.mkdir(parents=True, exist_ok=True)

    def breakpoint_log2(species, ab):
        bps = tasks["breakpoints_mg_per_L"]
        return math.log2(bps.get(species, {}).get(ab, bps["default"][ab]))

    def predict(train_ids, test_ids):
        tr = labels[labels.isolate_id.isin(train_ids)]
        med = tr.groupby(["species", "antibiotic"]).log2_mic_point.median()
        overall = tr.groupby("antibiotic").log2_mic_point.median()       # fallback for unseen species
        rows = []
        for iid in test_ids:
            sp = meta.at[iid, "species"]
            for ab in tasks["antibiotics"]:
                mu = med.get((sp, ab), overall[ab])
                rows.append({"isolate_id": iid, "antibiotic": ab, "log2_mic_pred": mu,
                             "p_resistant": float(mu > breakpoint_log2(sp, ab))})
        return pd.DataFrame(rows)

    for split_file in sorted((PUBLIC / "splits").glob("*.json")):
        split = json.loads(split_file.read_text())
        parts = split["partitions"]
        name = split["split"]
        pred = predict(parts["train"] + parts["val"], parts["test"])
        pred.to_csv(OUT / f"predictions_{name}.tsv", sep="\t", index=False)

    (OUT / "submission.json").write_text(json.dumps(
        {"model": "species_median", "team": "example external team",
         "description": "Per species x antibiotic median training MIC; ignores the genome"}, indent=2))
    print(f"submission written to {OUT}")


if __name__ == "__main__":
    main()
