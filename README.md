# AMRBench pilot: building the benchmark before the data exists

This folder is a working, end-to-end version of AMRBench that runs on **synthetic data shaped like the real deliverables**. When the UCL isolates are sequenced, you point the config at the real data, and the same code builds the splits, trains the three baselines (k-mer, pan-genome, GNN), scores them on T1–T7 and produces the leaderboard.

> All numbers produced here come from mock data. They test the *machinery*, not the science.

📖 **Full documentation: [the project wiki](wiki/Home.md)**. It covers the data contract, the submission format, metric definitions, the decision log, the roadmap and more.

```bash
cd /Users/k2259885/GIGABIOME_TECH/AMRBench/pilot
.venv/bin/amrbench run-all          # ~2 minutes: mock data -> leaderboard
```

---

## The big picture

A benchmark is three things: **a fixed dataset with fixed splits**, **fixed task definitions and metrics**, and **a harness that scores anyone's predictions the same way**. The code is organised as a pipeline where each stage only talks to the next through files:

```
  ┌──────────────┐   data contract   ┌───────────┐  release/public  ┌──────────┐  submissions/  ┌─────────┐
  │ sequencing + │ ───────────────▶  │   build   │ ───────────────▶ │  models  │ ─────────────▶ │ harness │──▶ leaderboard
  │ QC pipeline  │  (files, fixed    │ labels,   │  (labels for     │ (ours or │  (prediction   │ (reads  │
  │  OR mock     │   formats)        │ splits)   │  train/val only) │ external)│   files only)  │ private │
  └──────────────┘                   └───────────┘                  └──────────┘                │ labels) │
                                          │ release/private (test labels) ─────────────────────▶└─────────┘
```

The two seams that make the scaffold reusable:
1. **The data contract** ([amrbench/data.py](amrbench/data.py)): mock and real data enter through the same files. Swapping data means changing what *writes* those files, not what reads them.
2. **The submission format** ([amrbench/harness.py](amrbench/harness.py)): the harness never imports a model. It scores TSV files, so Phylax, a genomic LLM or an external team's model all go through the same door.

---

## Step-by-step

### Step 0: Setup (done already)
A virtual environment lives in `.venv/` (numpy, pandas, scikit-learn, PyTorch). The package is installed in editable mode, so `amrbench <step>` works. The tests run with `.venv/bin/python -m pytest tests`.

### Step 1: The data contract and the mock data: `amrbench generate-mock`
**Concept.** First decide *exactly* what the upstream pipeline must deliver. Everything downstream depends on it.

| File | Real-world producer | Content |
|---|---|---|
| `metadata.tsv` | UCL + QC report | isolate_id, species, collection_site, date, sequencing_batch, coverage, N50, contamination, qc_pass |
| `mic.tsv` | UCL broth microdilution | isolate_id, antibiotic, mic **as the lab reports it** (`0.25`, `<=0.0625`, `>64`) |
| `assemblies/<id>.fasta` | Unicycler | contigs |
| `pangenome/gene_presence_absence.Rtab` | Panaroo | gene family × isolate 0/1 |
| `pangenome/gene_calls.tsv` | Bakta GFF + Panaroo IDs | isolate, contig, start, end, strand, gene_family |
| `annotation/amrfinder_calls.tsv` | AMRFinderPlus | known determinants found per isolate (used to score T7) |
| `reference/known_determinants.tsv` | CARD / AMRFinder DB | which gene/mutation affects which antibiotic |

**Why the mock is more than plumbing.** [amrbench/mock/generate.py](amrbench/mock/generate.py) simulates a small world with a **known answer**:
- acquired resistance genes on plasmids shared across species, gyrA/parC point mutations, and species-level intrinsic resistance;
- **deliberate traps**: plasmid backbone genes (`repA`, `tra`, `mob`) that travel with resistance genes but cause nothing, and clonal lineages that share both marker genes and resistance;
- lab noise (±1 dilution), censored MICs, missing tests and QC failures.

Only the mock has a `truth/` folder. That's what lets us check that the benchmark measures what we think it measures.

### Step 2: Validate the inputs: `amrbench validate`
This checks every file against the contract (required columns, every MIC is a two-fold dilution, every QC-passed isolate has an assembly and a pan-genome entry, unknown species/antibiotics). **This is the first command to run on real data.** It will catch formatting mismatches long before any model runs.

### Step 3: Build the benchmark: `amrbench build` ([amrbench/build.py](amrbench/build.py))

**3a. Labels.** MICs are converted to log2 dilution steps (`0.25 mg/L → -2`), because "within ±1 dilution" is then just `|pred - true| ≤ 1`. The R/S call uses the breakpoint (R if MIC > breakpoint). Censored values need care: `<=8` with an 8 mg/L breakpoint is surely S, but `<=16` is genuinely unknown and stays unlabelled for T5.

**3b. Clone clusters.** Each genome is sketched as a set of 21-mers (like Mash). Pairwise Jaccard distances are computed, and isolates within a species are clustered (threshold `clone_distance_threshold`). On the mock data this recovers all 40 planted lineages exactly.

**3c. Leak-free splits.** `StratifiedGroupKFold` keeps every clone cluster entirely in train, val or test, while balancing species × collection site. This produces:
- `main` = 70% train / 10% val / 20% test, used for T1–T5 and T7;
- `T6_<species>` = the main split with one species removed from training. Its test set is that species' isolates *from the main test set*, so the in-distribution and held-out scores are compared on the very same isolates and no test labels leak.

**3d. Audit.** Every split manifest (`release/public/splits/*.json`) carries a leakage audit. Here is what it looked like on the mock data:

| split | clone clusters shared train↔test | test genomes with a near-identical training genome |
|---|---|---|
| clone-grouped (`main`) | 0 | 0 / 114 |
| naive random (internal) | 34 | **109 / 114** |

**3e. Public vs private.** `release/public/` holds metadata, train/val labels, split manifests and task definitions; this is what you would publish. `release/private/labels_all.tsv` holds the test labels, which **only the harness reads**.

### Step 4: Features: `amrbench features` ([amrbench/features.py](amrbench/features.py))
Each baseline sees the genome differently. The features are computed once, without labels:

| Model | Input | Can see SNPs? | Can see gene gain/loss? | Can see genomic context? |
|---|---|---|---|---|
| k-mer | hashed canonical 15-mers of the assembly | yes (one SNP = 15 new k-mers) | yes | no |
| pan-genome | Panaroo presence/absence + species | **no** | yes | no |
| GNN | gene graph (node = gene + allele, edge = neighbouring genes) | yes (via allele) | yes | yes |

### Step 5: Baselines: `amrbench train` ([amrbench/models/](amrbench/models/))
All three implement one interface ([base.py](amrbench/models/base.py)): `fit(train, val)`, `predict(ids)` and `explain(ids)`.
- **k-mer**: one ridge regression per antibiotic (solved in its dual form, because genomes ≪ k-mers).
- **pan-genome**: one random forest per antibiotic.
- **GNN**: plain-PyTorch message passing. All species and antibiotics are trained jointly, with early stopping on validation.

**Calibration (needed for ECE in T5–T7).** Each model predicts a log2 MIC. The validation set tells us how wrong it typically is (σ per antibiotic), which turns the prediction into a probability: `P(R) = P(true MIC > breakpoint)`. This is what the validation partition is for: tuning and calibration, never final scores.

Each model writes a **submission**: `submissions/<model>/predictions_main.tsv`, `attributions_main.tsv` and `predictions_T6_<species>.tsv`.

### Step 6: The harness: `amrbench evaluate` ([amrbench/harness.py](amrbench/harness.py), [amrbench/metrics.py](amrbench/metrics.py))
The harness reads a submission folder, validates it (columns, no duplicates, only test isolates, probabilities in [0,1]) and scores every task defined in the config:

| Task | Metric |
|---|---|
| T1–T4 | essential agreement (within ±1 dilution, handling censored MICs) with a bootstrap 95% CI over isolates; T4 per species |
| T5 | F1, balanced accuracy, sensitivity/specificity, plus the clinical **very major error** rate (called S, actually R) and **major error** rate |
| T6 | EA drop (in-distribution − held-out) and ECE on the held-out species |
| T7 | top-k recall of known determinants (per isolate × antibiotic), split by mechanism, plus ECE |

**How an external team takes part:** [examples/external_submission_species_median.py](examples/external_submission_species_median.py) uses only `release/public/`, writes prediction files, and is scored with `amrbench evaluate --submission submissions/species_median`. It is also the **floor baseline** (it ignores the genome entirely).

### Step 7: Leaderboard: `amrbench leaderboard`
Results go to `results/leaderboard.md`, and the per-task details go to `results/<model>/scores.json`.

---

## What the pilot already taught us (mock data, but the effects are real)

| model | T1 EA % | T3 EA % | T5 F1 | T6 drop pp | T7 top-5 recall | T7 recall, **mutations** | top-5 slots taken by plasmid decoys |
|---|---|---|---|---|---|---|---|
| GNN | 68.8 | 82.3 | 0.91 | 15.1 | 0.86 | 0.06 | 19% |
| k-mer | 75.3 | 94.1 | 0.96 | 9.5 | 0.86 | 0.03 | 28% |
| pan-genome | 80.5 | 100 | 0.90 | 18.4 | 0.85 | 0.00 | 41% |
| species median (floor) | 41.6 | 44.1 | 0.68 | n/a | n/a | n/a | n/a |

1. **Leaky splits inflate scores.** The same pan-genome model gets 94% EA on a random split but 83% on the clone-grouped split (`amrbench leakage-demo`). The 11-point gap is the over-optimism the proposal's R&D Challenge 2 is about, now measured.
2. **T6 behaves as intended.** *K. pneumoniae* (resistance on shared plasmids) transfers with roughly 0 pp drop. *E. faecium* (intrinsic resistance nobody else has) drops about 30 pp for every model. The leaderboard averages the two, so read the per-species values in `scores.json`.
3. **T7 headline recall hides a real problem.** Every model finds plasmid resistance genes (qnr: 40/40), but almost none credit the gyrA/parC mutations. For the k-mer model I checked that it does detect the mutation: the 15 new k-mers carry top-1% weights. Because the mutation is fixed within particular clones, though, clone marker genes get as much credit. The pan-genome model is also structurally blind to SNPs. On top of that, 19–41% of top-5 attributions are non-causal plasmid backbone genes. This is exactly the correlation-vs-causation gap that causaLens's work is meant to close, and the mock gives a known answer to test it against. **Recommendation: report T7 per mechanism, not just one recall number.**
4. **Small test sets mean wide confidence intervals.** T1 has 13 test isolates in the pilot, and its CI spans about ±10 points. With 1,500 real isolates there will be about 40 test isolates per species, still small enough that leaderboard differences of a few points may not be significant. Always report the CIs.

---

## Switching to real data: checklist

1. **Make the pipeline emit the contract.** The QC pipeline (`AMRBench_QC`: Unicycler, Bakta, CheckM2) already covers assemblies and annotation. It still needs **Panaroo** (→ `gene_presence_absence.Rtab`, plus `gene_data.csv` → `gene_calls.tsv`) and **AMRFinderPlus** (→ `amrfinder_calls.tsv`, including POINT mutations). A small adapter script will convert their native outputs into the contract files.
2. **The MIC export from UCL** → long format `isolate_id, antibiotic, mic`, keeping the lab's `<=`/`>` notation.
3. **Edit [config/benchmark.yaml](config/benchmark.yaml):** `data_root`, the real ≥15-antibiotic panel, the **EUCAST breakpoint table** (the current values are placeholders), and the tested MIC range.
4. Run `amrbench validate`, fix whatever it reports, then run `build features train evaluate leaderboard` (skip `generate-mock`).
5. **Re-tune `clone_distance_threshold`** by comparing clusters to MLST sequence types (add ST to the metadata).
6. **Scale the k-mer step.** Real genomes are about 1,000× larger than the mock genomes, and 2^18 hash bins would saturate. Use a dedicated k-mer/unitig tool (e.g. unitig-caller, or KMC filtered to variable k-mers) and write the same sparse matrix. The pan-genome model and the GNN scale as they are (the GNN would benefit from a GPU and PyTorch Geometric).
7. **Add Phylax** (and any genomic-LLM or multimodal baseline) by having it write a submission folder. No harness changes are needed.

## Open design decisions to agree with the team
- **Censored MICs**: the pilot trains on a point estimate (`>64` → 128). Censored-likelihood (interval) regression is the principled option.
- **EUCAST "I" (susceptible, increased exposure) and ATU**: the pilot uses only the R> breakpoint.
- **Site leakage**: the pilot *stratifies* by site. A stricter variant would hold out a whole site.
- **T7 resolution**: gene-level (pilot) or variant-level (e.g. gyrA S83L). Expert-panel review sits outside the harness.
- **Test-set size and repeated splits**, given the CI widths above.

## Code map
```
config/benchmark.yaml        everything that defines the benchmark (tasks, breakpoints, splits, models)
amrbench/data.py             data contract: readers + validation
amrbench/mock/generate.py    synthetic dataset with planted ground truth
amrbench/build.py            labels, clone clusters, splits, audit, public/private release
amrbench/features.py         k-mers, pan-genome matrix, gene graphs
amrbench/models/             base API + kmer.py, pangenome.py, gnn.py
amrbench/baselines.py        trains baselines -> submission folders
amrbench/harness.py          scores submission files against hidden labels
amrbench/metrics.py          EA, F1/VME, ECE, top-k recall, bootstrap CI
amrbench/leaderboard.py      results -> leaderboard.md
amrbench/experiments.py      internal leakage demo
examples/                    external-team submission example (floor baseline)
tests/                       unit tests for parsing, labels, metrics, k-mers, ridge
```
