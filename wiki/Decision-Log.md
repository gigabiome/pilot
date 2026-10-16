# Decision log

This page records benchmark design decisions so that partners (Gigabiome, causaLens, UCL) and later readers know *what* was decided and *why*. Add a row when something changes, and don't delete superseded rows: mark them **superseded**.

## Decided

| ID | Decision | Rationale | Status |
|---|---|---|---|
| D01 | The benchmark models are a **GNN, a k-mer model and a pan-genome model**, plus Phylax as the named baseline | Three complementary views of the genome (sequence, gene content, gene context). Replaces the application's list of model classes. | decided by the project lead |
| D02 | The harness scores **files**, never imports models | Any language or framework can take part; our baselines get no special treatment (R&D Challenge 3) | adopted in pilot |
| D03 | MICs are handled on the **log2 dilution scale**; EA is ±1 dilution with explicit censoring rules | Matches the Appendix Q9 metric; censored lab values are scored fairly | adopted in pilot |
| D04 | R/S rule: **R if MIC > R-breakpoint** (EUCAST "R >") | Simple and standard; the "I" category is still open (O2) | adopted in pilot |
| D05 | **Clone-grouped splits**: k-mer Jaccard distance, average-linkage clusters within species, StratifiedGroupKFold stratified on species × site | Prevents clonal leakage (R&D Challenge 2); measured at 11 pp inflation otherwise | adopted in pilot; threshold to re-tune (O4) |
| D06 | **T6 splits are derived from the main split** (held-out species' test = its main-test isolates) | No test-label leakage; in- and out-of-distribution scores on identical isolates | adopted in pilot |
| D07 | **Test labels are private** and read only by the harness | Benchmark integrity, as stated in the application | adopted in pilot |
| D08 | **Missing predictions count as wrong**; coverage is reported | Prevents cherry-picking easy pairs | adopted in pilot |
| D09 | **Bootstrap 95% CIs over isolates** for EA | Pairs within an isolate are correlated | adopted in pilot |
| D10 | T5 also reports **VME and ME**; T7 also reports **recall per mechanism** | Clinical relevance; the pilot showed the aggregate recall hides mutation failures | adopted in pilot |
| D11 | A **species-median floor** is included on the leaderboard | Separates genomic signal from species identity | adopted in pilot |
| D12 | Baselines calibrate P(R) from **validation-set residuals** | Gives every regression model a probability for ECE; participants may use their own method | adopted in pilot |

## Open

| ID | Question | Options | Notes |
|---|---|---|---|
| O1 | How to train on and score **censored MICs** | point estimate (current: `>v` → v+1) · censored/interval likelihood · exclude | affects all MIC tasks |
| O2 | **EUCAST "I"** (susceptible, increased exposure), ATU, and the **EUCAST version** | R vs non-R · three-class · report both | pin the version in config |
| O3 | **Site leakage** | stratify by site (current) · hold out a whole site as an extra split | only 2 sites, which limits hold-out |
| O4 | **Clone threshold for real data** | tune k-mer Jaccard against MLST · use PopPUNK · use Mash | check cluster ↔ ST agreement |
| O5 | **T7 resolution and expert review** | gene-level (current) · variant-level (e.g. gyrA S83L) · protocol for the expert panel | variant level needs SNP features in all models |
| O6 | **Test-set size / repeated splits** | single 20% split · larger test share · repeated grouped CV for the baselines | per-species CIs about ±10 pp in the pilot |
| O7 | **Freezing splits across sequencing batches** | build once at M9 · assign new batches without moving old isolates · one final build at M12 | must be settled before any interim results |
| O8 | **External data policy** for submissions | allowed if declared and test-overlap removed (draft in [Submission format](Submission-Format.md)) · not allowed | needs public-genome overlap checks |
| O9 | **k-mer tooling at scale** | unitig-caller · KMC + filter · larger hash space | see [Real-data migration](Real-Data-Migration.md) |
| O10 | **Which species to hold out in T6** | one per mechanism type (e.g. one with shared plasmid resistance, one with intrinsic resistance) · rotate all 7 | the pilot shows the choice drives the result |
| O11 | **Which species × drug pairs are scored** | all · exclude EUCAST "expected resistant phenotypes" and pairs without breakpoints · report them separately | intrinsic pairs inflate the scores |
| O12 | **Gene-family harmonisation across species** | Panaroo per species · Panaroo across all species · cross-species clustering (MMseqs2) | determines whether pan-genome/GNN features transfer in T6 |
