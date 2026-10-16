# Tasks and metrics

The tasks are defined in `config.tasks` and scored in [`harness.py`](../amrbench/harness.py). The metric implementations are in [`metrics.py`](../amrbench/metrics.py) and are unit-tested.

All tasks except T6 are scored from **one** prediction file on the main test split. The harness filters it by species or label as each task requires.

## Tasks

| ID | Kind | Scored pairs | Primary metric | Target |
|---|---|---|---|---|
| T1 | `mic` | test × *S. aureus* | essential agreement | ≥ 0.85 |
| T2 | `mic` | test × *A. baumannii* | essential agreement | ≥ 0.88 |
| T3 | `mic` | test × *P. aeruginosa* | essential agreement | ≥ 0.75 |
| T4 | `mic`, `per_species` | test × *Enterobacter*, *Enterococcus*, *E. coli*, *Klebsiella* | EA per species; the leaderboard shows the minimum | ≥ 0.75 each |
| T5 | `binary` | all test pairs with a known R/S | F1, sensitivity | F1 ≥ 0.88, sensitivity ≥ 0.90 |
| T6 | `cross_species` | test isolates of each held-out species | EA drop (pp), ECE | drop ≤ 15, ECE ≤ 0.10 |
| T7 | `attribution` | resistant test pairs that carry a known determinant | top-5 recall, ECE | recall ≥ 0.80, ECE ≤ 0.10 |

To add a task, add an entry to `config.tasks` with an existing `kind`. A new `kind` needs a scorer function in `harness.py`.

## Metrics

### Essential agreement (EA): T1–T4, T6
The prediction is rounded to the nearest dilution, r = round(pred). A pair is correct if:

| censor | Correct when |
|---|---|
| `=` | \|r − v\| ≤ 1 |
| `<=` (true ≤ v) | r ≤ v + 1 |
| `>` (true ≥ v + 1) | r ≥ v |

**A missing prediction counts as wrong**, and the `coverage` field reports how many were missing.

**Confidence interval:** a 95% bootstrap CI, resampling whole **isolates** (500 resamples). Pairs from the same isolate are correlated, so resampling pairs would understate the uncertainty.

### Binary metrics: T5
R is the positive class. The R/S call is `p_resistant ≥ 0.5` if probabilities are supplied for every pair; otherwise it is `log2_mic_pred > breakpoint`.

| Metric | Definition |
|---|---|
| sensitivity | TP / (TP + FN) |
| specificity | TN / (TN + FP) |
| balanced accuracy | mean of sensitivity and specificity |
| F1 | 2·precision·sensitivity / (precision + sensitivity) |
| **very major error (VME)** | FN / (TP + FN): called S, actually R. The clinically dangerous error. |
| **major error (ME)** | FP / (TN + FP): called R, actually S |

VME and ME follow the AST device-evaluation convention (ISO 20776-2). `per_antibiotic_f1` is also reported.

### Expected calibration error (ECE): T5, T6, T7
Predictions are binned by `p_resistant` into 10 equal-width bins:

ECE = Σ_bins (n_bin / n) · |observed R rate − mean predicted P(R)|

0 means perfectly calibrated. ECE is only computed when every pair has a `p_resistant`.

### Cross-species drop: T6
For each held-out species S:
- `ea_in_distribution` = EA of the **main** predictions on S's test isolates (the model saw S in training);
- `ea_held_out` = EA of the `T6_S` predictions on the same isolates (the model never saw S);
- `ea_drop_pp` = 100 × (in − held out).

The leaderboard shows the mean drop and mean ECE across held-out species. **Always read the per-species values**, because they can differ greatly; see [Pilot findings](Pilot-Findings.md).

### Top-k recall: T7
For every test pair that is (a) resistant and (b) carries at least one known determinant for that antibiotic (according to `amrfinder_calls.tsv` ∩ `known_determinants.tsv`):

recall@k = |top-k attributed features ∩ determinants present| / |determinants present|

The reported value is the mean over pairs (k = `top_k`, 5 by default). Intrinsic resistance is not scored because there is no gene to point at.

Also reported:
- `top5_recall_by_mechanism`: the fraction of determinants found, split by `acquired_gene` and `point_mutation`;
- `top5_decoy_rate_MOCK_ONLY`: the fraction of top-k slots taken by known non-causal plasmid backbone genes. This only exists for mock data, because it needs `truth/`.

The **expert-panel agreement** part of T7 is a manual review outside the harness.

## `scores.json` structure
```
{
  "submission": {...from submission.json...},
  "benchmark": "AMRBench-pilot", "version": "0.1.0",
  "tasks": {
    "T1": {"name", "target", "essential_agreement", "ci95", "n_pairs", "n_isolates", "coverage", "per_antibiotic"},
    "T4": {..., "per_species": {species: {...}}, "min_species_ea"},
    "T5": {"n", "balanced_accuracy", "f1", "sensitivity", "specificity",
           "very_major_error_rate", "major_error_rate", "ece", "per_antibiotic_f1"},
    "T6": {species: {"ea_in_distribution", "ea_held_out", "ea_drop_pp", "ece_held_out", "n_pairs"}},
    "T7": {"ece", "top5_recall", "n_pairs", "top5_recall_by_mechanism", "top5_decoy_rate_MOCK_ONLY"}
  }
}
```
`results/<model>/test_pairs_scored.tsv` has every test pair joined with its prediction and EA outcome, for error analysis.
