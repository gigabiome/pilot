# Pilot findings

> These results come from **mock data** (seed 42, 543 QC-passed isolates). They show how the benchmark *behaves*, not how real models will perform. The effects below come from mechanisms that also exist in real clinical data.

## Leaderboard (mock)

| Model | T1 EA % | T2 EA % | T3 EA % | T4 min EA % | T5 F1 | T5 sens | T5 VME % | T6 drop pp | T6 ECE | T7 top-5 recall | T7 recall, mutations | T7 ECE | decoy rate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GNN | 68.8 | 71.4 | 82.3 | 70.4 | 0.906 | 0.895 | 10.5 | 15.1 | 0.047 | 0.859 | 0.064 | 0.041 | 0.19 |
| k-mer | 75.3 | 86.7 | 94.1 | 70.4 | 0.960 | 0.953 | 4.7 | 9.5 | 0.104 | 0.855 | 0.032 | 0.045 | 0.28 |
| pan-genome | 80.5 | 87.8 | 100 | 75.9 | 0.896 | 0.949 | 5.1 | 18.4 | 0.096 | 0.849 | 0.000 | 0.052 | 0.41 |
| species median (floor) | 41.6 | 57.1 | 44.1 | 28.3 | 0.678 | 0.652 | 34.8 | 10.4 | 0.124 | — | — | 0.260 | — |
| **Target (App. Q9)** | 85 | 88 | 75 | 75 | 0.88 | 0.90 | | ≤ 15 | ≤ 0.10 | 0.80 | | ≤ 0.10 | |

To reproduce: `amrbench run-all`, then run the floor example (see [Getting started](Getting-Started.md)).

## Finding 1: Leaky splits inflate scores
The same pan-genome model gets **94.3% EA on a random split and 82.8% on the clone-grouped split** (`amrbench leakage-demo`). In the random split, 109 of 114 test genomes had a near-identical training genome; in the clone-grouped split, none did.

**Implication:** without clone-grouped splits, a published benchmark would overstate performance by about 11 points here. This quantifies R&D Challenge 2.

## Finding 2: The genome adds a lot, and the floor is needed to see it
Every genome model beats the species-median floor on T1–T4, by 14–56 points (usually 30 or more). The floor's T5 very-major-error rate is 35%. Without the floor, we could not tell real genomic signal apart from "knowing the species".

## Finding 3: T6 behaves as intended, but the mean hides it
| Held-out species | k-mer drop | pan-genome drop | GNN drop |
|---|---|---|---|
| *K. pneumoniae* (resistance on shared plasmids) | −8.3 | 6.7 | −1.7 |
| *E. faecium* (intrinsic resistance unique to the species) | 27.3 | 30.1 | 31.8 |

Resistance carried on transferable elements generalises across species; species-specific intrinsic resistance cannot. The leaderboard mean (9.5–18.4 pp) mixes the two.

**Implication:** report T6 per held-out species, and choose the real held-out species deliberately ([Decision log](Decision-Log.md) O10).

## Finding 4: T7 headline recall hides a causal-attribution failure
Overall top-5 recall is about 0.85 for all models, which is above the 0.80 target. Broken down by mechanism:
- **acquired genes: recall 1.00** for every model (e.g. qnr found 40/40);
- **point mutations (gyrA/parC): 0.00–0.06.**

Why this happens:
- **The pan-genome model** cannot see SNPs at all.
- **The k-mer model** *does* detect the gyrA mutation: the 15 mutation-specific k-mers sit in the top 1% of coefficients. But in *E. coli* the mutation is present in 25/25, 10/11 and 9/9 isolates of lineages L1–L3 and absent from L0, L4 and L5. A correlative model therefore gives clone-marker genes as much credit as the causal gene.
- **Decoys:** 19–41% of top-5 slots go to plasmid backbone genes (`repA`, `tra`, `mob`) that travel with resistance genes but cause nothing.

**Implications:**
- T7 must be reported per mechanism. This is now done in the harness.
- This is exactly the gap causal methods (causaLens) should close. The mock data provides a known answer to test them against before any expert-panel time is spent.

## Finding 5: Uncertainty is large at per-species test sizes
T1 has 13 test isolates in the pilot, and its 95% CIs span about ±10 points (e.g. pan-genome 80.5% [71.4, 88.5]). With 1,500 real isolates, there will be about 40 test isolates per species.

**Implication:** leaderboard differences of a few points per species may not be significant. Always show the CIs, and consider repeated splits for the baselines ([Decision log](Decision-Log.md) O6).

## Finding 6: Calibration via the validation set works
Every baseline reaches T5/T7 ECE ≤ 0.06 with the simple validation-residual calibration. On the unseen species *E. faecium*, the k-mer and pan-genome models degrade to ECE 0.14–0.15, above the 0.10 target, while the GNN stays at 0.03. Calibration under distribution shift therefore separates the models, which is exactly what T6 should expose.

## What the mock cannot tell us
- Real model rankings. The mock is simpler and much smaller than real genomes.
- Real lab noise, batch effects, or MIC distributions per drug.
- Whether 0.12 is the right clone threshold for real data.
- How well AMRFinderPlus annotation matches the real causal truth.
