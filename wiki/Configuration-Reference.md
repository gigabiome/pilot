# Configuration reference

Everything that defines the benchmark lives in [`config/benchmark.yaml`](../config/benchmark.yaml). Relative paths are resolved from the project root. Use another file with `--config`.

## Top level
| Key | Pilot value | Meaning |
|---|---|---|
| `benchmark_name`, `version` | `AMRBench-pilot`, `0.1.0` | written into manifests and scores. Bump `version` whenever the splits or labels change. |
| `seed` | 42 | drives the mock data, the splits and the models |
| `data_root` | `data/mock` | **the input folder ([data contract](Data-Contract.md)); point it at the real data** |
| `work_dir` | `work` | intermediate files (features, distances) |
| `release_dir` | `release` | `public/` + `private/` |
| `submissions_dir`, `results_dir` | `submissions`, `results` | |
| `species` | the 7 ESKAPE-E species | exact names, as in `metadata.tsv` |
| `antibiotics` | 6 drugs | **replace with the real ≥15-drug panel** |

## `mic_panel`
| Key | Pilot value | Meaning |
|---|---|---|
| `min_log2`, `max_log2` | −4, 6 | tested range (0.0625–64 mg/L); the mock generator censors outside it |

## `breakpoints`
```yaml
breakpoints:
  default: {ciprofloxacin: 0.5, ...}       # mg/L, R if MIC > value
  Staphylococcus aureus: {ciprofloxacin: 2}  # species overrides
```
⚠️ These are **placeholders**. Replace them with EUCAST values. See [Labels and breakpoints](Labels-and-Breakpoints.md).

## `qc`
| Key | Value | Meaning |
|---|---|---|
| `require_qc_pass` | true | exclude isolates with `qc_pass = False` |

## `splits`
| Key | Pilot value | Meaning |
|---|---|---|
| `test_fraction` | 0.2 | approximate test share (whole clusters are kept together) |
| `val_fraction` | 0.1 | approximate validation share |
| `sketch_k`, `sketch_bits` | 21, 20 | k-mer size and hash width for the genomic-distance sketch |
| `clone_distance_threshold` | 0.12 | Jaccard distance used to cut clone clusters. **Re-tune it on real data.** |
| `stratify_by` | `[species, collection_site]` | metadata columns kept balanced across partitions |
| `cross_species_holdout` | *K. pneumoniae*, *E. faecium* | one T6 split per entry |

## `tasks`
One entry per task: `name`, `kind` (`mic` / `binary` / `cross_species` / `attribution`), `target`, and kind-specific keys:
- `mic`: `species` (list) and optionally `per_species: true`;
- `attribution`: `top_k`.

See [Tasks and metrics](Tasks-and-Metrics.md).

## `models`
| Model | Keys |
|---|---|
| `kmer` | `k` (15), `n_bits` (18), `min_prevalence` (0.02), `alphas` |
| `pangenome` | `n_estimators` (300), `min_samples_leaf` (2), `max_features` (0.3), `min_prevalence` (0.02) |
| `gnn` | `hidden` (64), `layers` (3), `epochs` (150), `lr` (0.003), `batch_size` (32), `patience` (25) |

## `mock`
This section is only used by `generate-mock`. See [Mock data](Mock-Data.md#knobs).
`isolates_per_species`, `n_lineages`, `n_core_genes`, `n_accessory_genes`, `gene_length`, `mic_noise_sd`, `missing_rate`, `qc_fail_rate`, `sites`, `intrinsic`.
