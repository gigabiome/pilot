# Splits and leakage

**Leakage** happens when information about the test set reaches training. In clinical AMR data the main source is **clonality**: outbreaks yield near-identical genomes. If siblings sit on both sides of the train/test line, a model can memorise the clone and look much better than it will on a new patient. Preventing this is R&D Challenge 2. The code is in [`build.py`](../amrbench/build.py).

## Step 1: Genomic distance (Mash-like)
Each assembly is turned into its set of canonical 21-mers (`sketch_k`), hashed into 2^20 bins (`sketch_bits`). The distance is **Jaccard**: 1 − |A∩B| / |A∪B|. It is computed for all pairs with one sparse matrix product. The result is saved to `work/distances.npy`.

## Step 2: Clone clusters
Within each species, average-linkage hierarchical clustering cuts the tree at `clone_distance_threshold` (0.12 in the pilot). Every isolate gets a `clone_cluster` ID such as `Escherichia_coli_C003`, which is written into `release/public/metadata.tsv`.

**Check on mock data:** 40 clusters, and every cluster corresponds to exactly one of the 40 planted lineages.

## Step 3: The main split
`StratifiedGroupKFold` takes `groups = clone_cluster`, so a cluster is never divided between partitions, and `y = species|collection_site`, so the partitions stay balanced.
- One fold of 5 becomes **test** (20%).
- From the rest, one fold becomes **val** (10% of the total).
- The remainder is **train** (~70%).

Splitting whole clusters makes the final fractions approximate: the pilot gives 375 / 54 / 114 isolates.

## Step 4: T6 cross-species splits
For each species in `cross_species_holdout`, the split is derived **from the main split**:
- train / val = main train / val **minus** that species;
- test = that species' isolates **in the main test set**.

This has two benefits. No test label of the main split is ever used for training, and the in-distribution model (main split) and the held-out model are scored on the **very same isolates**.

## Step 5: Audit
Every split gets an audit, which is stored in its manifest and printed by `amrbench build`:

| Field | Meaning | Healthy value |
|---|---|---|
| `clone_clusters_shared_train_test` | clusters present in both train and test | 0 |
| `test_isolates_with_near_identical_train_genome` | test isolates whose nearest same-species training genome is closer than the threshold | 0 |
| `median_nearest_train_distance` | typical novelty of the test genomes | well above the threshold |
| `species_counts`, `site_counts` | balance per partition | similar proportions |
| `resistance_rate` | % R per antibiotic per partition | similar rates |

Pilot results:

| Split | Shared clusters | Near-identical test genomes | Median nearest distance |
|---|---|---|---|
| `main` (clone-grouped) | 0 | 0 / 114 | 0.203 |
| naive random (internal) | 34 | **109 / 114** | 0.028 |

`amrbench leakage-demo` shows what this means in practice: the same pan-genome model scores **94.3%** essential agreement on the random split and **82.8%** on the clone-grouped split.

## Manifests (`release/public/splits/<name>.json`)
```json
{
  "benchmark": "AMRBench-pilot", "version": "0.1.0", "split": "main",
  "created": "2026-10-16T13:36:06+00:00",
  "method": "clone-grouped, stratified by species+collection_site",
  "seed": 42,
  "params": { "test_fraction": 0.2, "val_fraction": 0.1, "clone_distance_threshold": 0.12, "...": "..." },
  "labels_sha256": "c302464f…",
  "partitions": { "train": ["AMRB…"], "val": ["…"], "test": ["…"] },
  "audit": { "clone_clusters_shared_train_test": 0, "...": "..." }
}
```
`labels_sha256` ties the manifest to the exact label file it was built against.

## Rules for using the splits
- **Train** is for fitting; **val** is for tuning and calibration; **test** is used only by the harness.
- Never tune on test, even informally ("I tried 5 settings and kept the best test score").
- The random split exists only for the internal demo. It is never released or reported as a result.

## Open questions
Whether to hold out whole collection sites, and how to choose the threshold for real data, are in the [Decision log](Decision-Log.md).
