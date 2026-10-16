# Project overview

## The problem
AI models can predict antibiotic resistance from a genome in hours instead of the 24–48 h that culture-based testing takes. But there is no standard, open way to compare these models. Public resources (BV-BRC, NCBI AMR, CARD) are fragmented by species, rarely have clinical **MIC** values, and come with no evaluation infrastructure. Most public genomes with phenotypes carry only binary R/S tags, gathered under heterogeneous protocols.

## What AMRBench delivers
1. **A homogeneous clinical dataset**: ≥1,500 ESKAPE-E isolates from UCL microbiology and the Royal Free Hospital. All are tested by EUCAST broth microdilution against ≥15 antibiotics in one lab, and all are sequenced with Illumina short reads.
2. **An open benchmark package** (CC-BY 4.0): task definitions, leak-free splits released as JSON manifests, and baseline results.
3. **A pip-installable evaluation harness** (Apache 2.0) and a public leaderboard.
4. **A licensed curated dataset** for commercial AI training (tiered access).

## Consortium and work packages

| WP | Lead | Scope | Months |
|---|---|---|---|
| WP1 | R. Jalan (UCL) | Sample collection, EUCAST MIC testing, biobanking, MIxS metadata | M-3 – M12 |
| WP2 | F. Clasen (Gigabiome) | WGS (NovoGene UK, 4 batches), assembly, QC, Nextflow pipeline | M1 – M12 |
| WP3 | M. Sipos (causaLens) | Task design, splits, baselines, harness, external validation | M3 – M12 |
| WP4 | F. Clasen (Gigabiome) | Open release, licensing portal, manuscript, engagement | M9 – M12 |

**The pilot is WP3 work started early.** It lets the benchmark and harness be designed and tested before WP1/WP2 deliver data.

## The seven tasks (from Appendix Q9)

| ID | Task | Primary metric | Target |
|---|---|---|---|
| T1 | MIC prediction, *S. aureus* | essential agreement (±1 dilution) | ≥ 85% |
| T2 | MIC prediction, *A. baumannii* | essential agreement | ≥ 88% |
| T3 | MIC prediction, *P. aeruginosa* | essential agreement | ≥ 75% |
| T4 | MIC prediction, *Enterobacter*, *Enterococcus*, *E. coli*, *Klebsiella* | essential agreement per species | ≥ 75% each |
| T5 | Binary R/S at the EUCAST breakpoint | F1, sensitivity, balanced accuracy | F1 ≥ 0.88, sensitivity ≥ 90% |
| T6 | Cross-species generalisation | EA drop vs in-distribution, ECE | drop ≤ 15 pp, ECE ≤ 0.10 |
| T7 | Causal feature attribution & calibration | top-k recall of known determinants, expert agreement | top-5 recall ≥ 0.80, expert agreement ≥ 80% |

The full definitions are on [Tasks and metrics](Tasks-and-Metrics.md).

## What has changed since the application
The application allows tasks and named baselines to be adapted. The changes so far:

| In the application | Now | Notes |
|---|---|---|
| Baselines: genomic LLM, GNN, multimodal clinical-genomic model, Phylax | **Three benchmark models: GNN, k-mer model, pan-genome model** (plus Phylax as the named baseline) | All three are implemented in the pilot. Any further model (e.g. a genomic LLM) plugs in through the [submission format](Submission-Format.md). |
| — | A **species-median floor baseline** | Added so every result can be compared against "ignoring the genome". |
| T7 reported as one recall number | T7 recall is also reported **per mechanism** (acquired gene vs point mutation) | Prompted by a pilot finding; see [Pilot findings](Pilot-Findings.md). |
| T5 metrics: F1, balanced accuracy, sensitivity/specificity | **Very major error** and **major error** rates added | These are the standard clinical AST error rates. |

Record further changes in the [Decision log](Decision-Log.md).

## Three R&D challenges, and where the pilot addresses them
1. **Causal feature attribution.** T7 plus planted confounders in the mock data, which give a known answer to test attribution against. See [Mock data](Mock-Data.md).
2. **Leak-free stratified splits.** Clone-grouped, stratified splits with a published leakage audit. See [Splits and leakage](Splits-and-Leakage.md).
3. **Harness usability across heterogeneous models.** A file-based submission format and the external-team example. See [Submission format](Submission-Format.md).
