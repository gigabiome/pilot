# Submission format

This is the public contract between AMRBench and **any** model: ours, Phylax, or an external team's. The harness ([`harness.py`](../amrbench/harness.py)) reads only these files. You can use any language, framework or hardware.

## What a participant receives
Everything in `release/public/`:

| File | Use |
|---|---|
| `metadata.tsv` | species, site and clone cluster of every isolate |
| `labels_trainval.tsv` | MIC labels for train + val isolates ([columns](Labels-and-Breakpoints.md#label-file-columns)) |
| `splits/main.json` | train / val / test isolate IDs |
| `splits/T6_<species>.json` | cross-species splits |
| `tasks.json` | tasks, antibiotics, breakpoints, MIC panel |

They also get the genome inputs described in the [data contract](Data-Contract.md): assemblies, pan-genome and annotations. They never receive test labels.

## What a participant returns
```
submissions/<name>/
├── submission.json                 required
├── predictions_main.tsv            required: T1–T5, T7
├── attributions_main.tsv           optional: T7 recall
└── predictions_T6_<species>.tsv    optional: T6, one per held-out species
```

### `submission.json`
```json
{"model": "my_model", "team": "My Lab", "description": "one sentence on the method"}
```

### `predictions_main.tsv` (and `predictions_T6_*.tsv`)
One row per test isolate × antibiotic:

| Column | Required | Notes |
|---|---|---|
| `isolate_id` | ✅ | must be in that split's test set |
| `antibiotic` | ✅ | as in `tasks.json` |
| `log2_mic_pred` | ✅* | predicted log2 MIC (mg/L). *Alternatively supply `mic_pred` in mg/L and it is log2-converted. |
| `p_resistant` | recommended | P(MIC > breakpoint), in [0, 1]. Needed for ECE; without it T5 uses `log2_mic_pred` vs the breakpoint. |

```
isolate_id  antibiotic     log2_mic_pred  p_resistant
AMRB00003   ciprofloxacin  -2.241         0.1185
```
Each T6 file must come from a model trained **only** on that split's train/val partitions.

### `attributions_main.tsv`
For T7, give ranked genomic features per isolate × antibiotic. Provide at least `top_k` rows (5).

| Column | Notes |
|---|---|
| `isolate_id`, `antibiotic` | |
| `rank` | 1 = most important |
| `feature` | gene family name, matching the pan-genome / `known_determinants.tsv` vocabulary |
| `score` | optional |

## Validation rules (a submission is rejected with `SubmissionError` if…)
- a required column is missing;
- there is a duplicate `isolate_id` × `antibiotic` row;
- an isolate is not in that split's test set;
- `p_resistant` falls outside [0, 1].

Rows for unknown antibiotics are ignored. **Missing pairs are scored as wrong**, and coverage is reported.

## Scoring a submission
```bash
.venv/bin/amrbench evaluate --submission submissions/<name>
```
```bash
.venv/bin/amrbench leaderboard
```

## Rules of the benchmark
1. Train only on train; tune and calibrate only on val.
2. Don't use test labels in any form, including other people's results on the test set.
3. *(Draft, pending [Decision log](Decision-Log.md) O8)* External data (e.g. public BV-BRC genomes) is allowed if declared in `description`. Isolates that overlap AMRBench test isolates must be removed.
4. The T6 models must never see the held-out species.

## Worked example
[`examples/external_submission_species_median.py`](../examples/external_submission_species_median.py) reads only `release/public/` and writes a complete submission for the main and T6 splits. It doesn't submit attributions, which shows that the harness handles a missing optional file.
