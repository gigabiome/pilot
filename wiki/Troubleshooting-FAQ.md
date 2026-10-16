# Troubleshooting & FAQ

## Errors

| Message | Cause | Fix |
|---|---|---|
| `amrbench: command not found` | the venv is not activated | use `.venv/bin/amrbench`, or `source .venv/bin/activate` |
| `ContractError: metadata.tsv is missing required columns` | the input doesn't match the [data contract](Data-Contract.md) | rename or add the columns; extra columns are fine |
| `ContractError: MIC '3' is not a two-fold dilution step` | a non-dilution value in `mic.tsv` | check the lab export (units? typos?) |
| `ContractError: N isolates have no assembly` | a FASTA is missing or misnamed | the file must be `assemblies/<isolate_id>.fasta` |
| `ContractError: … missing from the pan-genome matrix` | Panaroo was run on a different isolate set | re-run Panaroo on all QC-passed isolates |
| `SubmissionError: … not in this test split` | predictions for train/val isolates, or the wrong split | predict exactly the `test` IDs of that manifest |
| `SubmissionError: duplicate isolate_id x antibiotic rows` | duplicate rows | one row per pair |
| `FileNotFoundError: … work/features/ids.txt` | `train` was run before `features` | run the steps in order, or `run-all` |
| `KeyError` in `BenchData.y` | a model requested labels for unknown isolates | use only the IDs from the split manifests |

## FAQ

**Why is there a `truth/` folder? Won't models cheat?**
It exists only for mock data, and model code never reads it. The harness uses it for one mock-only diagnostic (decoy rate). Real data won't have it.

**Why do my numbers differ from the wiki?**
Check that you haven't changed `seed` or anything in `config/benchmark.yaml`. Everything is deterministic. If the config is unchanged, a clean run (`rm -rf data/mock work release submissions results`, then `run-all`) reproduces the [Pilot findings](Pilot-Findings.md).

**The split sizes aren't exactly 70/10/20. Is that a bug?**
No. Whole clone clusters are assigned to partitions, so the fractions are approximate. See [Splits and leakage](Splits-and-Leakage.md).

**Why is T6 drop sometimes negative?**
The held-out model did slightly better than the in-distribution model on that species. With small test sets (10 *Klebsiella* isolates in the pilot), this is noise.

**Can I tune on the test set if I don't report it?**
No. Any decision informed by test scores leaks information. Use the validation partition.

**Why ridge, random forest and a small GNN? Aren't there better models?**
Baselines should be simple, transparent and reproducible, so that they set a clear reference. Stronger models (Phylax, genomic LLMs) are expected to beat them. See [Baseline models](Baseline-Models.md#adding-a-model).

**Can I submit without attributions or T6 files?**
Yes. T7 recall and T6 then show "not submitted"; everything else is scored.

**Where do I change which species are held out, or the target numbers?**
In [`config/benchmark.yaml`](../config/benchmark.yaml). See [Configuration reference](Configuration-Reference.md).

**How do I edit this wiki?**
Edit the Markdown files in `pilot/wiki/`. New pages need a link in `Home.md` and `_Sidebar.md`.
