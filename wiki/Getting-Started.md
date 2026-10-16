# Getting started

## Requirements
- Python ≥ 3.10 (developed on 3.13)
- About 50 MB of disk for the mock run
- No GPU needed for the pilot

## Install
The project already has a virtual environment at `pilot/.venv`. On a fresh machine:

```bash
cd /Users/k2259885/GIGABIOME_TECH/AMRBench/pilot
python3 -m venv .venv
.venv/bin/pip install -e ".[gnn,dev]"
```

This installs the `amrbench` command and the dependencies: numpy, pandas, scipy, scikit-learn, PyYAML, PyTorch and pytest.

## Run everything
```bash
.venv/bin/amrbench run-all
```
This takes about 2 minutes on a laptop. The final lines print the leaderboard and the leakage demo.

## Run step by step (recommended the first time)
Each step reads the previous step's files and writes its own, so you can stop after any step and inspect the output.

| # | Command | Reads | Writes | Look at |
|---|---|---|---|---|
| 1 | `amrbench generate-mock` | config | `data/mock/` | `data/mock/metadata.tsv`, `mic.tsv`, `truth/ground_truth.json` |
| 2 | `amrbench validate` | `data/mock/` | — (prints a report) | the console |
| 3 | `amrbench build` | `data/mock/` | `release/public/`, `release/private/`, `work/distances.npy` | `release/public/splits/main.json` → `audit` |
| 4 | `amrbench features` | `data/mock/`, `release/public/metadata.tsv` | `work/features/` | — |
| 5 | `amrbench train` | `release/public/`, `work/features/` | `submissions/<model>/` | `predictions_main.tsv`, `attributions_main.tsv` |
| 6 | `amrbench evaluate` | `submissions/`, `release/private/` | `results/<model>/scores.json` | `scores.json` |
| 7 | `amrbench leaderboard` | `results/*/scores.json` | `results/leaderboard.md` | `leaderboard.md` |
| 8 | `amrbench leakage-demo` | `release/private/` | `results/leakage_demo.tsv` | the console |

Useful options:
```bash
.venv/bin/amrbench train --models kmer gnn
```
```bash
.venv/bin/amrbench evaluate --submission submissions/species_median
```
```bash
.venv/bin/amrbench build --config path/to/other.yaml
```

## Try the external-team path
```bash
.venv/bin/python examples/external_submission_species_median.py
```
```bash
.venv/bin/amrbench evaluate --submission submissions/species_median
```
```bash
.venv/bin/amrbench leaderboard
```

## Run the tests
```bash
.venv/bin/python -m pytest tests
```
The 16 tests cover MIC parsing, censored R/S calls, every metric, k-mer strand-independence, the SNP → k-mer effect, Jaccard distance, and the dual ridge solver (checked against scikit-learn).

## Generated folders
`data/mock/`, `work/`, `release/`, `submissions/` and `results/` are all regenerated from the config and code. They are listed in `.gitignore`, and deleting them is always safe.

Next: [Architecture](Architecture.md).
