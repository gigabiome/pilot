# Data contract

The data contract is the set of files the upstream pipeline (WP1 + WP2) must deliver into `data_root`. It is defined and validated in [`amrbench/data.py`](../amrbench/data.py). **The mock generator and the real pipeline produce the same files.** Everything downstream depends only on this page.

```
<data_root>/
├── metadata.tsv
├── mic.tsv
├── assemblies/<isolate_id>.fasta
├── pangenome/gene_presence_absence.Rtab
├── pangenome/gene_calls.tsv
├── annotation/amrfinder_calls.tsv
└── reference/known_determinants.tsv
```
All tables are tab-separated with a header row.

---

## `metadata.tsv`: one row per isolate
**Producer:** UCL metadata manifest (MIxS) joined with the WP2 QC report.

| Column | Type | Example | Notes |
|---|---|---|---|
| `isolate_id` | string, unique | `AMRB00001` | primary key everywhere |
| `species` | string | `Escherichia coli` | must match a name in `config.species` exactly |
| `collection_site` | string | `Royal Free` | used to stratify the splits |
| `collection_date` | ISO date | `2027-01-30` | |
| `sequencing_batch` | string | `B3` | NovoGene batch (B1–B4) |
| `coverage` | float | `82.6` | mean read depth |
| `n50` | int | `228531` | assembly N50 (bp) |
| `contamination_pct` | float | `0.62` | e.g. from CheckM2 |
| `qc_pass` | bool | `True` | accepts `True/False`, `1/0`, `yes/no`. Failing isolates are excluded when `qc.require_qc_pass` is set. |

Extra columns are allowed and ignored. MLST sequence type (`st`) is recommended, to sanity-check the clone clusters.

## `mic.tsv`: one row per measured isolate × antibiotic
**Producer:** UCL broth-microdilution results.

| Column | Example | Notes |
|---|---|---|
| `isolate_id` | `AMRB00001` | must exist in `metadata.tsv` |
| `antibiotic` | `meropenem` | lower case, must match `config.antibiotics` |
| `mic` | `0.25`, `<=0.0625`, `>64` | **exactly as the lab reports it**, in mg/L |

- Accepted prefixes are `<=`, `<`, `>`, `>=`, `≤` and `≥`. The parser normalises them; see [Labels and breakpoints](Labels-and-Breakpoints.md).
- Values must be two-fold dilution steps (…, 0.25, 0.5, 1, 2, …). Anything else raises `ContractError`.
- Untested pairs are simply absent (no row). Missing tests are normal.

## `assemblies/<isolate_id>.fasta`
**Producer:** Unicycler. Standard FASTA, one record per contig. The contig name is the first word of the header and must match `gene_calls.tsv`.

## `pangenome/gene_presence_absence.Rtab`
**Producer:** Panaroo (the Roary format is identical).

```
Gene          AMRB00001  AMRB00002  AMRB00003
aac6_mock     1          1          0
gyrA          1          1          1
```
Rows are gene families, columns are isolates, and values are 0/1. Every QC-passed isolate must have a column.

## `pangenome/gene_calls.tsv`: one row per predicted gene
**Producer:** Bakta GFF coordinates, mapped to Panaroo family IDs (from Panaroo's `gene_data.csv`).

| Column | Example | Notes |
|---|---|---|
| `isolate_id` | `AMRB00001` | |
| `contig` | `AMRB00001_contig1` | must match the FASTA header |
| `start`, `end` | `341`, `640` | 1-based, inclusive |
| `strand` | `+` / `-` | |
| `gene_family` | `gyrA` | same identifier as the Rtab row |

This file is used to build the [gene graphs](Features.md) and to map k-mers back to genes for T7.

## `annotation/amrfinder_calls.tsv`: known determinants found per isolate
**Producer:** AMRFinderPlus (run with `--organism` so that point mutations are reported).

| Column | Example | Notes |
|---|---|---|
| `isolate_id` | `AMRB00001` | |
| `gene_family` | `blaESBL_mock`, `gyrA` | must use the same names as `known_determinants.tsv` |
| `element_subtype` | `AMR` / `POINT` | |

This is the reference for scoring T7: which known determinants each isolate actually carries.

## `reference/known_determinants.tsv`: which feature affects which drug
**Producer:** curated from the CARD / AMRFinder databases, then reviewed by the clinical microbiologists.

| Column | Example |
|---|---|
| `antibiotic` | `ciprofloxacin` |
| `gene_family` | `gyrA` |
| `mechanism` | `acquired_gene` / `point_mutation` |

One determinant can affect several antibiotics; it gets one row per antibiotic.

---

## Validation (`amrbench validate`)
**Fatal problems** (raise `ContractError`):
- missing required columns, or duplicate `isolate_id`;
- MIC rows for isolates with no metadata;
- an MIC that isn't a two-fold dilution;
- a QC-passed isolate without an assembly or a pan-genome column.

**Reported, not fatal:**
- species or antibiotics not listed in the config (ignored);
- grid coverage, censoring counts and measurements per antibiotic.

## Mock-only folder
`truth/` (`ground_truth.json`, `lineages.tsv`) exists only for mock data. Nothing in the benchmark requires it. The harness uses it for one diagnostic, if present: the decoy rate in T7.
