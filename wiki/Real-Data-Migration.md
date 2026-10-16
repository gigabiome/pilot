# Real-data migration

The goal is that switching from mock to real data changes **only the inputs and the config**. This page lists what has to happen, in order.

## Checklist

### A. Upstream: make the pipeline produce the data contract
- [ ] **Assemblies.** Unicycler output → `assemblies/<isolate_id>.fasta`. This already exists in `AMRBench_QC`.
- [ ] **QC metadata.** Merge the UCL metadata with QC metrics (coverage, N50, CheckM2 contamination, pass/fail) → `metadata.tsv`. Agree the `qc_pass` rule (proposal: ≥30× coverage, N50 > 50 kb, plus a contamination limit).
- [ ] **Pan-genome.** Add **Panaroo** to the pipeline, using Bakta GFFs as input. `gene_presence_absence.Rtab` is used as is.
- [ ] **Gene calls.** Write an adapter: Bakta GFF coordinates + Panaroo `gene_data.csv` (gene → cluster) → `gene_calls.tsv`.
- [ ] **AMR calls.** Add **AMRFinderPlus** with `--organism`, so point mutations are reported. Write an adapter → `amrfinder_calls.tsv`, normalising gene symbols to the gene-family vocabulary.
- [ ] **Known determinants.** Build `known_determinants.tsv` from the AMRFinderPlus/CARD databases for the chosen antibiotics, and have the clinical microbiologists review it.
- [ ] **MIC.** UCL exports a long table `isolate_id, antibiotic, mic`, keeping `<=`/`>` exactly as reported.

> The adapters belong in a new `amrbench/adapters/` module. Each is a small script that converts a native tool output into a contract file.

### B. Configuration (`benchmark.yaml`)
- [ ] `data_root: data/real`
- [ ] `antibiotics`: the real ≥15-drug panel, with names matching `mic.tsv`.
- [ ] `breakpoints`: **the EUCAST table**, species-specific, with the EUCAST version pinned and recorded in `version`.
- [ ] `mic_panel`: the real tested ranges (these may differ per drug; extend the config if so).
- [ ] `cross_species_holdout`: chosen deliberately ([Decision log](Decision-Log.md) O10).
- [ ] Bump `version`.

### C. Run and check
```bash
.venv/bin/amrbench validate
```
Fix everything it reports, then:
```bash
.venv/bin/amrbench build
```
- [ ] **Re-tune `clone_distance_threshold`.** Add MLST STs to the metadata and check that clusters roughly match STs, or sub-split large STs during outbreaks. As an alternative, consider PopPUNK clusters.
- [ ] Read the audit: 0 shared clusters, balanced species/site counts, and similar resistance rates per partition.

### D. Scale the features
| Component | Mock | Real (~5 Mb genomes, ~1,500 isolates) | Action |
|---|---|---|---|
| k-mers | 45 kb, 2^18 bins | ~5M k-mers per genome; 2^18 would saturate | Use **unitig-caller** or KMC, keep only variable k-mers/unitigs (2–98% prevalence), and write the same sparse matrix plus gene mapping |
| distance sketch | 2^20 bins | OK, or use **Mash** directly | optional |
| pan-genome RF | ~400 genes | ~10–30k families | fine; maybe pre-filter by prevalence |
| GNN | ~130 nodes per graph | ~5,000 nodes per graph | GPU (A10G), PyTorch Geometric, mini-batching |
| dual ridge | n = 375 | n ≈ 1,050 | fine (an n × n solve) |

### E. Baselines and release
- [ ] `amrbench features train evaluate leaderboard`
- [ ] Phylax writes a submission folder (no harness changes).
- [ ] Invite external teams with `release/public/` + inputs + [Submission format](Submission-Format.md). This is the M9 gate.
- [ ] Before publishing: anonymisation check of the metadata (UK GDPR / NHS IG), the licence files (CC-BY 4.0 data, Apache 2.0 harness), and Zenodo DOIs.

## Dress rehearsal (recommended before batch 1)
Run the pipeline on a **public dataset with real genomes and MICs** (e.g. a BV-BRC subset of one or two ESKAPE-E species) to shake out the adapters and the scaling issues before UCL data arrives. Only the adapters and the config change; the benchmark code stays the same.

## Incremental batches
Sequencing arrives in four batches (M3, M6, M9, M12). Running `build` on a growing dataset **reshuffles** the splits each time. Decide how the splits are frozen ([Decision log](Decision-Log.md) O7) before results are reported on any interim split.
