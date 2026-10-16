"""Turn genomes into model inputs. Three "views" of the same genome:

  k-mers      : every length-k substring of the assembly (alignment-free, sees SNPs)
  pan-genome  : which gene families are present (sees gene gain/loss, blind to SNPs)
  gene graph  : genes as nodes, genomic neighbours as edges, with an allele label per
                node (sees genes, their alleles and their genomic context)

Features are computed WITHOUT labels, so they can be computed once for all isolates.
"""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from amrbench.data import assembly_path, load_gene_calls, load_gene_presence, read_fasta

_CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
_GOLDEN = np.uint64(0x9E3779B97F4A7C15)


# ----------------------------------------------------------------------------- k-mers

def kmer_bins(seq: str, k: int, n_bits: int) -> np.ndarray:
    """Hashed bins of all CANONICAL k-mers in `seq` (unique, sorted).

    Canonical = min(k-mer, reverse complement), so a gene scores the same whichever
    strand the assembler put it on. Hashing into 2^n_bits bins keeps the feature space
    fixed-size (real tools: KMC/Jellyfish + a fixed vocabulary, or unitigs/DBGWAS).
    k <= 31 so a k-mer fits in 64 bits (2 bits per base).
    """
    codes = _CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    if len(codes) < k:
        return np.empty(0, dtype=np.int64)
    win = np.lib.stride_tricks.sliding_window_view(codes, k)
    valid = (win != 255).all(axis=1)                      # skip windows containing N
    shifts = (2 * np.arange(k - 1, -1, -1)).astype(np.uint64)
    w = win.astype(np.uint64)
    fwd = (w << shifts).sum(axis=1, dtype=np.uint64)
    rc = ((np.uint64(3) - w[:, ::-1]) << shifts).sum(axis=1, dtype=np.uint64)
    canon = np.minimum(fwd, rc)[valid]
    with np.errstate(over="ignore"):
        h = (canon * _GOLDEN) >> np.uint64(64 - n_bits)
    return np.unique(h.astype(np.int64))


def kmer_matrix(ids: list[str], root: Path, k: int, n_bits: int) -> sp.csr_matrix:
    """Binary isolate x hashed-kmer presence matrix (sparse)."""
    rows, cols = [], []
    for r, iid in enumerate(ids):
        bins = np.unique(np.concatenate([kmer_bins(s, k, n_bits) for s in read_fasta(assembly_path(root, iid)).values()]))
        rows.append(np.full(len(bins), r))
        cols.append(bins)
    data = np.ones(sum(len(c) for c in cols), dtype=np.float32)
    return sp.csr_matrix((data, (np.concatenate(rows), np.concatenate(cols))), shape=(len(ids), 2 ** n_bits))


def jaccard_distance(X: sp.csr_matrix) -> np.ndarray:
    """Pairwise Jaccard distance between rows of a binary sparse matrix (Mash-like)."""
    X = (X > 0).astype(np.float32)
    inter = (X @ X.T).toarray()
    size = np.asarray(X.sum(axis=1)).ravel()
    union = size[:, None] + size[None, :] - inter
    return 1.0 - inter / np.maximum(union, 1)


def gene_sequences(root: Path, iid: str, calls: pd.DataFrame) -> list[tuple[str, str, str, int]]:
    """[(gene_family, dna_sequence, contig, start)] for one isolate (calls = its rows of gene_calls.tsv)."""
    contigs = read_fasta(assembly_path(root, iid))
    sub = calls.sort_values(["contig", "start"])
    return [(r.gene_family, contigs[r.contig][r.start - 1:r.end], r.contig, r.start) for r in sub.itertuples()]


# ----------------------------------------------------------------------------- gene graphs

def gene_graph(genes: list[tuple[str, str, str, int]]) -> dict:
    """One graph per genome: nodes = genes, edges = adjacent genes on the same contig.

    Each node carries its gene family and an allele hash (identical sequence -> same
    allele). With real data the allele label could be replaced by a protein language
    model embedding (e.g. ESM) of the gene.
    """
    families = [g[0] for g in genes]
    alleles = [hashlib.md5(g[1].encode()).hexdigest()[:12] for g in genes]
    src, dst = [], []
    for i in range(len(genes) - 1):
        if genes[i][2] == genes[i + 1][2]:                 # same contig
            src += [i, i + 1]
            dst += [i + 1, i]
    return {"family": families, "allele": alleles, "edges": np.array([src, dst], dtype=np.int64)}


# ----------------------------------------------------------------------------- build all

def build_features(cfg, ids: list[str]) -> Path:
    root, out = cfg.path("data_root"), cfg.path("work_dir") / "features"
    out.mkdir(parents=True, exist_ok=True)
    km = cfg["models"]["kmer"]

    X = kmer_matrix(ids, root, km["k"], km["n_bits"])
    sp.save_npz(out / "kmer.npz", X)

    pa = load_gene_presence(root).reindex(ids).fillna(0).astype(np.uint8)
    pa.to_pickle(out / "pangenome.pkl")

    calls = dict(tuple(load_gene_calls(root).groupby("isolate_id")))
    graphs, gene_bins = {}, {}
    for iid in ids:
        genes = gene_sequences(root, iid, calls[iid])
        graphs[iid] = gene_graph(genes)
        # Per-gene k-mer bins: lets the k-mer model map important k-mers back to genes (T7).
        gene_bins[iid] = [(fam, kmer_bins(seq, km["k"], km["n_bits"])) for fam, seq, _, _ in genes]
    with open(out / "graphs.pkl", "wb") as fh:
        pickle.dump(graphs, fh)
    with open(out / "gene_kmer_bins.pkl", "wb") as fh:
        pickle.dump(gene_bins, fh)
    pd.Series(ids).to_csv(out / "ids.txt", index=False, header=False)
    return out


class FeatureStore:
    """Read-only access to the precomputed features, row-aligned to `ids`."""

    def __init__(self, cfg):
        d = cfg.path("work_dir") / "features"
        self.ids = pd.read_csv(d / "ids.txt", header=None)[0].tolist()
        self.row = {iid: i for i, iid in enumerate(self.ids)}
        self._dir = d
        self._cache: dict = {}

    def _load(self, name, fn):
        if name not in self._cache:
            self._cache[name] = fn()
        return self._cache[name]

    def kmer(self, ids: list[str]) -> sp.csr_matrix:
        X = self._load("kmer", lambda: sp.load_npz(self._dir / "kmer.npz").tocsr())
        return X[[self.row[i] for i in ids]]

    def pangenome(self, ids: list[str]) -> pd.DataFrame:
        return self._load("pa", lambda: pd.read_pickle(self._dir / "pangenome.pkl")).loc[ids]

    def graphs(self, ids: list[str]) -> list[dict]:
        g = self._load("graphs", lambda: pickle.load(open(self._dir / "graphs.pkl", "rb")))
        return [g[i] for i in ids]

    def gene_kmer_bins(self, iid: str):
        return self._load("gkb", lambda: pickle.load(open(self._dir / "gene_kmer_bins.pkl", "rb")))[iid]
