"""Baseline 3 - graph neural network on per-genome gene graphs.

Input : one graph per genome. Node = gene (family embedding + allele embedding),
        edge = neighbouring genes on the same contig.
Model : message passing (each layer: node <- f(node, mean of neighbours)), then
        mean+max pooling to a genome vector, then one output per antibiotic.
        All species and antibiotics are trained jointly (multi-task, cross-species).
Sees  : genes, alleles (so point mutations) and genomic context (e.g. which genes
        share a plasmid contig).
Explain: gradient x input on each node, aggregated per gene family.

Written in plain PyTorch so every step is visible. For the real data you would
likely switch to PyTorch Geometric and richer node features (e.g. ESM protein
embeddings) - the BaselineModel interface stays the same.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
from torch import nn

from amrbench.models.base import BaselineModel

UNK = 0


class GeneGraphNet(nn.Module):
    def __init__(self, n_families, n_alleles, n_out, hidden, layers):
        super().__init__()
        self.fam = nn.Embedding(n_families, hidden)
        self.allele = nn.Embedding(n_alleles, hidden)
        self.self_lin = nn.ModuleList(nn.Linear(hidden, hidden) for _ in range(layers))
        self.nbr_lin = nn.ModuleList(nn.Linear(hidden, hidden, bias=False) for _ in range(layers))
        self.head = nn.Sequential(nn.Linear(2 * hidden, hidden), nn.ReLU(), nn.Dropout(0.2), nn.Linear(hidden, n_out))

    def embed(self, fam, allele):
        return self.fam(fam) + self.allele(allele)

    def forward_from_embedding(self, h, edges, batch, n_graphs):
        src, dst = edges
        deg = torch.zeros(h.shape[0]).index_add_(0, dst, torch.ones_like(dst, dtype=torch.float)).clamp(min=1)
        for lin_s, lin_n in zip(self.self_lin, self.nbr_lin):
            msg = torch.zeros_like(h).index_add_(0, dst, h[src]) / deg[:, None]     # mean over neighbours
            h = torch.relu(lin_s(h) + lin_n(msg)) + h                              # residual
        counts = torch.zeros(n_graphs).index_add_(0, batch, torch.ones_like(batch, dtype=torch.float))
        mean = torch.zeros(n_graphs, h.shape[1]).index_add_(0, batch, h) / counts[:, None]
        mx = torch.full((n_graphs, h.shape[1]), -1e9).scatter_reduce(0, batch[:, None].expand_as(h), h, "amax")
        return self.head(torch.cat([mean, mx], dim=1))

    def forward(self, fam, allele, edges, batch, n_graphs):
        return self.forward_from_embedding(self.embed(fam, allele), edges, batch, n_graphs)


class GNNModel(BaselineModel):
    name = "gene_graph_gnn"
    description = "Message-passing GNN on gene-adjacency graphs (multi-task, cross-species)"

    # ----------------------------------------------------------------- graph -> tensors
    def _encode(self, g):
        fam = torch.tensor([self.fam_vocab.get(f, UNK) for f in g["family"]])
        allele = torch.tensor([self.allele_vocab.get((f, a), UNK) for f, a in zip(g["family"], g["allele"])])
        return fam, allele, torch.from_numpy(g["edges"])

    def _batch(self, encoded):
        fams, alleles, edges, batch, offset = [], [], [], [], 0
        for i, (f, a, e) in enumerate(encoded):
            fams.append(f); alleles.append(a); edges.append(e + offset)
            batch.append(torch.full((len(f),), i, dtype=torch.long))
            offset += len(f)
        return torch.cat(fams), torch.cat(alleles), torch.cat(edges, dim=1), torch.cat(batch), len(encoded)

    # ----------------------------------------------------------------- training
    def _fit(self, train_ids, val_ids):
        p = self.params
        torch.manual_seed(0)
        train_graphs = self.data.features.graphs(train_ids)
        # Vocabularies come from TRAINING genomes only; unseen genes/alleles map to UNK.
        self.fam_vocab, self.allele_vocab = {}, {}
        for g in train_graphs:
            for f, a in zip(g["family"], g["allele"]):
                self.fam_vocab.setdefault(f, len(self.fam_vocab) + 1)
                self.allele_vocab.setdefault((f, a), len(self.allele_vocab) + 1)

        ytr = self.data.y(train_ids).to_numpy(np.float32)
        self.y_mean, self.y_std = np.nanmean(ytr, 0), np.nanstd(ytr, 0) + 1e-6
        enc_tr = [self._encode(g) for g in train_graphs]
        enc_va = [self._encode(g) for g in self.data.features.graphs(val_ids)]
        Ytr = torch.tensor((ytr - self.y_mean) / self.y_std)
        Yva = torch.tensor((self.data.y(val_ids).to_numpy(np.float32) - self.y_mean) / self.y_std)

        self.net = GeneGraphNet(len(self.fam_vocab) + 1, len(self.allele_vocab) + 1,
                                len(self.data.antibiotics), p["hidden"], p["layers"])
        opt = torch.optim.AdamW(self.net.parameters(), lr=p["lr"], weight_decay=1e-4)
        best, best_state, wait = np.inf, None, 0
        rng = np.random.default_rng(0)
        for epoch in range(p["epochs"]):
            self.net.train()
            order = rng.permutation(len(enc_tr))
            for start in range(0, len(enc_tr), p["batch_size"]):
                b = order[start:start + p["batch_size"]]
                loss = self._loss(self.net(*self._batch([enc_tr[i] for i in b])), Ytr[b])
                opt.zero_grad(); loss.backward(); opt.step()
            self.net.eval()
            with torch.no_grad():
                val_loss = float(self._loss(self.net(*self._batch(enc_va)), Yva))
            if val_loss < best - 1e-4:
                best, wait = val_loss, 0
                best_state = {k: v.clone() for k, v in self.net.state_dict().items()}
            else:
                wait += 1
                if wait >= p["patience"]:
                    break
        self.net.load_state_dict(best_state)
        self.net.eval()
        self.epochs_trained = epoch + 1

    @staticmethod
    def _loss(pred, y):
        mask = ~torch.isnan(y)                    # untested antibiotics contribute nothing
        return nn.functional.huber_loss(pred[mask], y[mask])

    # ----------------------------------------------------------------- inference
    def _predict_log2(self, ids):
        enc = [self._encode(g) for g in self.data.features.graphs(ids)]
        out = []
        with torch.no_grad():
            for s in range(0, len(enc), 128):
                out.append(self.net(*self._batch(enc[s:s + 128])).numpy())
        pred = np.vstack(out) * self.y_std + self.y_mean
        return pd.DataFrame(pred, index=ids, columns=self.data.antibiotics)

    def _gene_scores(self, iid, antibiotic):
        if getattr(self, "_explained_iid", None) != iid:   # one backward pass per isolate, reused across antibiotics
            g = self.data.features.graphs([iid])[0]
            fam, allele, edges = self._encode(g)
            h0 = self.net.embed(fam, allele).detach().requires_grad_(True)
            batch = torch.zeros(len(fam), dtype=torch.long)
            out = self.net.forward_from_embedding(h0, edges, batch, 1)[0]
            self._node_attr = {}
            for j, ab in enumerate(self.data.antibiotics):
                grad, = torch.autograd.grad(out[j], h0, retain_graph=True)
                self._node_attr[ab] = (grad * h0).sum(1).detach().numpy() * self.y_std[j]
            self._explained_iid, self._families = iid, g["family"]
        scores: dict[str, float] = {}
        for f, s in zip(self._families, self._node_attr[antibiotic]):
            scores[f] = max(scores.get(f, -np.inf), float(s))
        return scores
