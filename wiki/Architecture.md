# Architecture

## Pipeline
```
  sequencing + QC pipeline          amrbench build              models                  harness
  (or mock generator)
 ┌────────────────────┐  data    ┌──────────────────┐ public ┌────────────────┐ files ┌──────────────┐
 │ assemblies, MIC,   │ contract │ labels           │───────▶│ k-mer          │──────▶│ validate     │
 │ metadata, Panaroo, │─────────▶│ clone clusters   │        │ pan-genome     │       │ score T1–T7  │──▶ leaderboard
 │ AMRFinderPlus      │          │ splits + audit   │        │ GNN            │       │              │
 └────────────────────┘          │ public / private │        │ external teams │       └──────▲───────┘
                                 └────────┬─────────┘        └────────────────┘              │
                                          └──────── private test labels ─────────────────────┘
```

Each stage talks to the next **only through files**. That has three consequences:
- any stage can be re-run, inspected or replaced on its own;
- external teams need nothing but files;
- moving to the real data touches only the first box.

## The two seams

### 1. The data contract
[`amrbench/data.py`](../amrbench/data.py) defines what enters the benchmark and validates it. The mock generator and the real pipeline both have to produce these files. See [Data contract](Data-Contract.md).

### 2. The submission format
[`amrbench/harness.py`](../amrbench/harness.py) scores prediction **files** and never imports model code. Our own baselines are scored through the same path as everyone else's. See [Submission format](Submission-Format.md).

## The public / private boundary
```
release/
  public/                  <- would be published
    metadata.tsv             QC-passed isolates + clone_cluster
    labels_trainval.tsv      labels for main-split train + val only
    splits/main.json         partitions + leakage audit
    splits/T6_<species>.json
    tasks.json               tasks, antibiotics, breakpoints, MIC panel
  private/                 <- never published
    labels_all.tsv           includes test labels; read only by the harness
    random_split_INTERNAL.json   leaky split for the internal leakage demo
```
Model code loads data through `BenchData` ([`models/base.py`](../amrbench/models/base.py)), which reads only `release/public/`. The one deliberate exception is `experiments.py` (`_PrivateData`), used for internal method checks.

## Directory layout
```
pilot/
├── config/benchmark.yaml     everything that defines the benchmark
├── amrbench/
│   ├── cli.py                `amrbench <step>` entry point
│   ├── config.py             config loading, breakpoint lookup
│   ├── data.py               data contract: readers + validation
│   ├── mock/generate.py      synthetic data with planted ground truth
│   ├── build.py              labels, clone clusters, splits, audit, release
│   ├── features.py           k-mers, pan-genome, gene graphs, FeatureStore
│   ├── models/
│   │   ├── base.py           BenchData + BaselineModel API + calibration
│   │   ├── kmer.py           ridge on hashed k-mers
│   │   ├── pangenome.py      random forest on gene presence/absence
│   │   └── gnn.py            message-passing GNN on gene graphs
│   ├── baselines.py          train baselines -> submission folders
│   ├── harness.py            score submissions
│   ├── metrics.py            EA, F1/VME/ME, ECE, top-k recall, bootstrap CI
│   ├── leaderboard.py        results -> leaderboard
│   └── experiments.py        internal leakage demo
├── examples/                 external-team submission (floor baseline)
├── tests/                    unit tests
└── wiki/                     this documentation
```

## Design principles
- **Configuration over code.** Species, antibiotics, breakpoints, split fractions, task targets and model hyper-parameters all live in [`benchmark.yaml`](Configuration-Reference.md).
- **Determinism.** A single `seed` drives the mock data, the splits and the models. Two runs give identical numbers.
- **Provenance.** Every split manifest records the config, seed, method, creation time and a SHA-256 of the label file.
- **Fail loudly at the boundary.** Contract violations raise `ContractError` and malformed submissions raise `SubmissionError`. Nothing silently coerces bad data.
