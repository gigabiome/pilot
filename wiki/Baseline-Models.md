# Baseline models

There are three benchmark baselines, a k-mer model, a pan-genome model and a GNN, plus a genome-blind floor. They all implement one interface and are scored through the same [submission format](Submission-Format.md) as external models.

## The model API ([`models/base.py`](../amrbench/models/base.py))

```python
class BaselineModel:
    def fit(self, train_ids, val_ids) -> self          # public labels only
    def predict(self, ids) -> DataFrame                # isolate_id, antibiotic, log2_mic_pred, p_resistant
    def explain(self, ids, top_k) -> DataFrame         # isolate_id, antibiotic, rank, feature, score
```
A subclass implements three methods:

| Method | Returns |
|---|---|
| `_fit(train_ids, val_ids)` | — |
| `_predict_log2(ids)` | DataFrame, ids × antibiotics, of predicted log2 MIC |
| `_gene_scores(iid, antibiotic)` | dict: gene family → contribution to raising the MIC |

The base class handles calibration, the output format and ranking (only positive scores are ranked).

`BenchData` is everything a model may see: public metadata, public labels (`y(ids)` returns `log2_mic_point`, with NaN where untested), breakpoints and the `FeatureStore`.

### Calibration: from a MIC prediction to P(resistant)
After `_fit`, the base class measures each antibiotic's error on the **validation** set:

σ = max(RMSE_val, 0.3)

Treating the true log2 MIC as Normal(μ, σ²), with breakpoint b:

P(R) = P(true > b) = 1 − Φ((b + 0.5 − μ) / σ)

The +0.5 puts the threshold halfway between the last S dilution and the first R dilution. This gives every regression model a calibrated probability, which T5–T7 (ECE) need.

## k-mer model: `kmer_ridge` ([`models/kmer.py`](../amrbench/models/kmer.py))
- **Input:** hashed k-mer presence. k-mers present in < 2% or > 98% of training genomes are dropped.
- **Model:** one ridge regression per antibiotic on log2 MIC. α is chosen from {1, 10, 100, 1000} by validation MAE.
- **Solver:** the dual (kernel) form, w = Xcᵀ(XcXcᵀ + αI)⁻¹(y − ȳ). It is exact (it matches scikit-learn to 1e-6 in the tests) and fast when genomes ≪ k-mers. K = XXᵀ is computed once per fit.
- **Explain:** for each gene in the isolate, the sum of the coefficients of its k-mers (the maximum over copies of the same family).
- **Strength:** sees mutations. **Weakness:** the weight of a causal k-mer is shared with every k-mer correlated with it (clone markers).

## Pan-genome model: `pangenome_rf` ([`models/pangenome.py`](../amrbench/models/pangenome.py))
- **Input:** gene presence/absence (with the same 2–98% prevalence filter) plus species one-hot.
- **Model:** one random forest per antibiotic (300 trees, `min_samples_leaf` = 2, `max_features` = 0.3).
- **Explain:** impurity importance of each gene the isolate carries. It doesn't show direction; SHAP is the natural upgrade.
- **Strength:** captures gene interactions; robust. **Weakness:** structurally blind to point mutations.

## GNN: `gene_graph_gnn` ([`models/gnn.py`](../amrbench/models/gnn.py))
- **Input:** gene graphs. Node embedding = family embedding + allele embedding. The vocabularies are built from **training genomes only**, and unseen families/alleles map to UNK.
- **Architecture:** 3 message-passing layers with hidden size 64. Each layer computes `h ← ReLU(W_self·h + W_nbr·mean(neighbours)) + h`, with a residual connection. Mean + max pooling then feeds an MLP head with one output per antibiotic.
- **Training:** multi-task (all antibiotics) and cross-species (all species). Targets are standardised per antibiotic, and the Huber loss is masked so untested pairs contribute nothing. AdamW (lr 0.003, weight decay 1e-4), batch size 32, up to 150 epochs, early stopping on validation loss (patience 25).
- **Explain:** gradient × input on each node's embedding, the maximum per gene family, scaled back to log2 units.
- Written in plain PyTorch so every operation is visible. For the real data, consider PyTorch Geometric and richer node features.

## Floor: `species_median` ([`examples/`](../examples/external_submission_species_median.py))
Predicts the median training MIC for each species × antibiotic and ignores the genome. Any genome model must beat it clearly. It is also the worked example of an external submission.

## Adding a model
1. Create `amrbench/models/<name>.py` with a `BaselineModel` subclass that implements `_fit`, `_predict_log2` and `_gene_scores`, and sets `name` and `description`.
2. Register it in `MODELS` in [`baselines.py`](../amrbench/baselines.py) (e.g. `"llm": GenomicLLMModel`).
3. Add its hyper-parameters under `models.<key>` in `benchmark.yaml`.
4. Run `amrbench train --models <key>`, then `amrbench evaluate` and `amrbench leaderboard`.

If the model lives outside this codebase (Phylax, or an external team's model), skip all of this and write the [submission files](Submission-Format.md) directly.

## Runtime (pilot, laptop CPU)
| Model | Main split | Each T6 split |
|---|---|---|
| k-mer | ~5 s | ~3 s |
| pan-genome | ~20 s | ~1 s |
| GNN | ~20 s | ~13 s |
