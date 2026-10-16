# Roadmap

The project started on 1 September 2026 (M1 = September 2026). Dates below are derived from Appendix Q14.

## Grant milestones relevant to the benchmark

| Month | Date | Event |
|---|---|---|
| M2 | Oct 2026 | **Pilot complete (this repository)** |
| M3 | Nov 2026 | Sequencing batch 1 · **Gate 1**: sample QC + batch 1 pass |
| M6 | Feb 2027 | Batch 2 · **Gate 2**: annotated dataset + Nextflow pipeline released (D2.3) · external teams invited |
| M9 | May 2027 | Batch 3 · **D3.1** curated dataset with stratified splits · **Gate 3**: external validation passed |
| M10 | Jun 2027 | **D3.2** baseline results across all tasks · D1.2 MIC profiles |
| M11 | Jul 2027 | **D3.3** pip-installable harness · **D3.4** external validation from ≥3 teams |
| M12 | Aug 2027 | Batch 4 · **D4.1** public release (CC-BY 4.0 / Apache 2.0) · manuscript |

## Proposed benchmark work plan

### Now → M3: get ready for batch 1
- [ ] Agree the open decisions with causaLens, prioritising O1, O2, O7, O11 and O12 ([Decision log](Decision-Log.md)).
- [ ] Add Panaroo + AMRFinderPlus to `AMRBench_QC`; write the `amrbench/adapters/`.
- [ ] Get the antibiotic panel and EUCAST breakpoints from UCL; update the config.
- [ ] Do a dress rehearsal on public genomes (BV-BRC subset) to test the adapters and scaling.
- [ ] Put the pilot under version control and set up CI that runs the tests.

### M3 → M6: first real data
- [ ] `validate` + `build` on batch 1; tune the clone threshold against MLST.
- [ ] Scale the k-mer features (unitigs); move the GNN to GPU / PyTorch Geometric.
- [ ] Interim baselines on interim splits (internal only, never reported as results).
- [ ] Freeze the submission format v1.0 and write the external-team guide (based on [Submission format](Submission-Format.md)).

### M6 → M9: external validation (Gate 3)
- [ ] Onboard ≥3 external teams; collect usability feedback on the harness.
- [ ] Freeze the splits according to O7; publish the manifests with audits (D3.1).
- [ ] Phylax submission; T7 expert-panel protocol.

### M9 → M12: results and release
- [ ] Baseline results across T1–T7 with CIs (D3.2).
- [ ] Package the harness for PyPI, write the docs, and build the public leaderboard (D3.3).
- [ ] Anonymisation review, licences, Zenodo DOIs, release (D4.1).

## Known technical debt in the pilot
- No adapters for Panaroo, Bakta or AMRFinderPlus yet.
- k-mer hashing will not scale to real genomes as is.
- The pan-genome explainer is direction-agnostic (planned upgrade: SHAP).
- A single MIC panel range for all drugs.
- The pilot folder is not yet a git repository.
