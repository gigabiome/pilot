# Features

`amrbench features` ([`features.py`](../amrbench/features.py)) computes three views of every QC-passed genome. **No labels are used**, so the features are computed once for all isolates and shared by every split.

| View | Built from | Sees SNPs | Sees gene gain/loss | Sees genomic context | Used by |
|---|---|---|---|---|---|
| k-mers | assembly FASTA | ✅ | ✅ | ❌ | k-mer model |
| pan-genome | Panaroo Rtab | ❌ | ✅ | ❌ | pan-genome model |
| gene graph | gene calls + FASTA | ✅ (via allele) | ✅ | ✅ | GNN |

## k-mers
- Every substring of length **k = 15** in each contig. k-mers never span contig boundaries, and windows containing `N` are skipped.
- **Canonical:** each k-mer is replaced by min(k-mer, reverse complement), so a gene gives the same features on either strand. This is unit-tested.
- 2-bit encoded into a 64-bit integer, multiplied by a fixed odd constant, and the top **n_bits = 18** bits are kept. Hashing into a fixed 2^18 space keeps the matrix a constant width.
- The output is a binary sparse matrix of isolates × 2^18 (`work/features/kmer.npz`).
- **A single SNP changes exactly k k-mers.** This is unit-tested, and it is why the k-mer model can in principle see point mutations.

For T7, the k-mer bins of every *gene* in every isolate are also stored (`gene_kmer_bins.pkl`), so important k-mers can be mapped back to genes.

> **Scaling warning.** Mock genomes are about 45 kb; real ones are about 5 Mb. At 2^18 bins the real k-mer space would saturate. See [Real-data migration](Real-Data-Migration.md).

## Pan-genome matrix
The Panaroo `gene_presence_absence.Rtab`, transposed to isolates × gene families and aligned to the isolate order (`pangenome.pkl`). Presence/absence is blind to SNPs by construction: a mutated gyrA is still "gyrA present".

## Gene graphs
There is one graph per genome (`graphs.pkl`):
- **Nodes** are the genes from `gene_calls.tsv`, in contig/position order. Each node has:
  - `family`: the Panaroo gene family;
  - `allele`: the first 12 hex characters of the MD5 of its DNA sequence (identical sequence means identical allele).
- **Edges** connect neighbouring genes on the same contig, in both directions. Genes on a plasmid contig are therefore connected to each other but not to the chromosome.

With real data the allele label could be replaced by a protein language-model embedding (e.g. ESM) of each gene.

## Files in `work/features/`
| File | Content |
|---|---|
| `ids.txt` | isolate order (the row index of every matrix) |
| `kmer.npz` | sparse k-mer presence matrix |
| `pangenome.pkl` | gene presence DataFrame |
| `graphs.pkl` | dict isolate → {family, allele, edges} |
| `gene_kmer_bins.pkl` | dict isolate → [(family, k-mer bins)] |

Models read these through `FeatureStore`, which caches each file on first use.
