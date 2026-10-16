# Mock data

[`amrbench/mock/generate.py`](../amrbench/mock/generate.py) writes a synthetic dataset in the [data contract](Data-Contract.md) format. It has two purposes:

1. **Plumbing.** Every stage of the benchmark can run before real data exists.
2. **Validation of the method.** Because we write the biology, we know the right answer. We can check that the splits prevent leakage, that T6 detects non-transferable resistance, and that T7 can tell causal features from correlated ones.

> Gene names ending in `_mock` are invented. Effect sizes are illustrative. Nothing here is a biological claim.

## The simulated world (default config)

| Element | Default | Purpose |
|---|---|---|
| Species | the 7 ESKAPE-E species | task structure as in the real data |
| Isolates | 80 per species (560 in total) | |
| Genome | 100 core + 40 accessory genes, 300 bp each, 40 bp spacers | small enough to run in seconds |
| Lineages | 6 per species, with uneven sizes | clonal population structure |
| Contigs | chromosome split into 2–5 contigs, plus one contig per plasmid | assembly-like fragmentation |

### Causal determinants (the ground truth)

| Antibiotic | Determinant | Mechanism | Effect (log2 dilutions) |
|---|---|---|---|
| meropenem | `blaCARB_mock` | acquired gene | +7 |
| ceftazidime | `blaESBL_mock` | acquired gene | +6 |
| ceftazidime | `blaCARB_mock` | acquired gene | +4 |
| gentamicin | `aac6_mock` | acquired gene | +5 |
| tetracycline | `tetA_mock` | acquired gene | +5 |
| trimethoprim | `dfrA_mock` | acquired gene | +6 |
| ciprofloxacin | `gyrA` | point mutation (position 249) | +3 |
| ciprofloxacin | `parC` | point mutation | +2 |
| ciprofloxacin | `qnr_mock` | acquired gene | +2 |

Effects add up. A single gyrA mutation lands near the breakpoint, while gyrA + parC is clearly resistant, which mirrors stepwise fluoroquinolone resistance.

**Intrinsic resistance** (species-level, +7, with no gene to point at): *E. faecium* to ceftazidime and trimethoprim; *P. aeruginosa* to tetracycline and trimethoprim; *A. baumannii* to trimethoprim.

### Plasmids: where the confounders come from

| Plasmid | Backbone (non-causal) | Cargo (causal) |
|---|---|---|
| pCARB | `repA_pCARB`, `traX_pCARB` | `blaCARB_mock`, `aac6_mock` |
| pESBL | `repA_pESBL`, `mobA_pESBL` | `blaESBL_mock`, `qnr_mock`, `dfrA_mock` |
| pTET | `repA_pTET` | `tetA_mock` |
| pCRYPT | `repA_pCRYPT`, `hyp_pCRYPT` | none (a decoy plasmid) |

Plasmid sequences are identical across species, which simulates horizontal transfer. When a plasmid is present, each cargo gene is present with 85% probability. Backbone genes are therefore strongly correlated with resistance but cause none of it. These are the **decoys** that T7 should not rank highly.

### Population structure
Each lineage has:
- a fixed accessory-gene profile (each gene present with p = 0.4), with 3% gain/loss per isolate;
- a per-plasmid carriage rate drawn from {0.03, 0.4, 0.9}, so resistance is clonal;
- a 35% chance that the gyrA mutation is fixed (95% of members carry it, compared with 3% elsewhere), and a parC mutation in 60% of those lineages;
- 30 neutral lineage SNPs, plus 5 private SNPs per isolate;
- a collection-site bias.

### Phenotype and lab realism
- **log2 MIC** = species wild-type baseline (2–4 dilutions below the breakpoint) + intrinsic effect + the sum of determinant effects + N(0, 0.6²) lab noise.
- Values are rounded to a dilution and **censored** at the panel limits (`<=0.0625`, `>64`).
- 4% of tests are missing, and 3% of isolates fail QC (low coverage and N50).

## Outputs
These are all the [data contract](Data-Contract.md) files, plus `truth/`:
- `ground_truth.json`: the determinants with their effects, the decoy list and the intrinsic table;
- `lineages.tsv`: the true lineage of each isolate, used to check the clone clustering.

## Knobs
These are set in the `mock:` section of [`benchmark.yaml`](Configuration-Reference.md): `isolates_per_species`, `n_lineages`, `n_core_genes`, `n_accessory_genes`, `gene_length`, `mic_noise_sd`, `missing_rate`, `qc_fail_rate`, `sites` and `intrinsic`. The determinant table and the plasmids are constants at the top of `generate.py`.

**Useful experiments:**
- Raise `mic_noise_sd` to see how essential agreement degrades as the lab gets noisier.
- Lower `isolates_per_species` to about 40 to see confidence intervals widen.
- Remove the clonal confounding of the mutations: in `_make_isolate`, give gyrA/parC the same mutation probability in every lineage instead of fixing them in some. Then check whether T7 recall for mutations recovers.
- Make plasmid carriage independent of lineage (one rate for every lineage in `_make_lineages`) and watch the decoy rate change.
