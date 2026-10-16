"""Generate a synthetic dataset that has the SAME SHAPE as the real AMRBench deliverables.

Why mock data is useful beyond plumbing: because we write the biology ourselves,
we know the ground truth, so we can check that the benchmark measures what we
think it measures. The simulated world contains the main traps of real AMR data:

  * Acquired resistance genes on plasmids, shared across species (horizontal transfer).
  * Point mutations in chromosomal genes (gyrA / parC -> ciprofloxacin).
  * Intrinsic, species-level resistance (no gene to point at).
  * CONFOUNDERS: plasmid backbone genes (repA, tra, mob) travel with the resistance
    genes, and clonal lineages share both accessory genes and resistance. A model can
    score well by learning these correlates - which T7 (attribution) should expose.
  * Clonal population structure - which naive random splits leak (see leakage-demo).
  * Lab noise (+-1 dilution), censored MICs, missing tests and QC failures.

All gene names ending in `_mock` are invented. Effect sizes are illustrative only.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from amrbench.config import Config

DNA = np.frombuffer(b"ACGT", dtype=np.uint8)
DETERMINANT_POS = 249   # codon position mutated for resistance mutations (think gyrA S83)

# Mobile elements shared across species: backbone genes are NOT causal, cargo genes are.
PLASMIDS = {
    "pCARB":  {"backbone": ["repA_pCARB", "traX_pCARB"], "cargo": ["blaCARB_mock", "aac6_mock"]},
    "pESBL":  {"backbone": ["repA_pESBL", "mobA_pESBL"], "cargo": ["blaESBL_mock", "qnr_mock", "dfrA_mock"]},
    "pTET":   {"backbone": ["repA_pTET"], "cargo": ["tetA_mock"]},
    "pCRYPT": {"backbone": ["repA_pCRYPT", "hyp_pCRYPT"], "cargo": []},
}

# The causal ground truth: (antibiotic, gene_family, mechanism, effect in log2 dilutions).
DETERMINANTS = [
    ("meropenem",     "blaCARB_mock", "acquired_gene", 7.0),
    ("ceftazidime",   "blaCARB_mock", "acquired_gene", 4.0),
    ("ceftazidime",   "blaESBL_mock", "acquired_gene", 6.0),
    ("gentamicin",    "aac6_mock",    "acquired_gene", 5.0),
    ("tetracycline",  "tetA_mock",    "acquired_gene", 5.0),
    ("trimethoprim",  "dfrA_mock",    "acquired_gene", 6.0),
    ("ciprofloxacin", "qnr_mock",     "acquired_gene", 2.0),
    ("ciprofloxacin", "gyrA",         "point_mutation", 3.0),
    ("ciprofloxacin", "parC",         "point_mutation", 2.0),
]
INTRINSIC_EFFECT = 7.0


def species_code(species: str) -> str:
    genus, epithet = species.split()[:2]
    return (genus[0] + epithet[:2]).lower()


def random_dna(rng: np.random.Generator, n: int) -> str:
    return DNA[rng.integers(0, 4, n)].tobytes().decode()


def point_mutate(seq: str, pos: int, rng: np.random.Generator | None = None, new_base: str | None = None) -> str:
    if new_base is None:
        choices = [b for b in "ACGT" if b != seq[pos]]
        new_base = choices[rng.integers(0, 3)]
    return seq[:pos] + new_base + seq[pos + 1:]


@dataclass
class SpeciesPool:
    name: str
    code: str
    core: list[tuple[str, str]]            # ordered chromosome backbone (family, wild-type sequence)
    accessory: list[tuple[str, str]]
    accessory_slot: np.ndarray             # accessory gene j is inserted after core gene accessory_slot[j]
    spacers: list[str]                     # intergenic spacer after each core slot
    mutant_base: dict[str, str]            # resistance-mutation base for gyrA / parC
    baseline: dict[str, float]             # wild-type modal log2 MIC per antibiotic


@dataclass
class Lineage:
    accessory_profile: np.ndarray
    plasmid_rate: dict[str, float]
    gyrA_fixed: bool
    parC_fixed: bool
    neutral_snps: list[tuple[int, int]]
    site_bias: float
    weight: float = 1.0


@dataclass
class Isolate:
    isolate_id: str
    species: str
    lineage: int
    contigs: dict[str, list[tuple[str, str]]] = field(default_factory=dict)   # contig -> [(family, seq)]
    carried: set[str] = field(default_factory=set)        # known determinants carried (gene families)
    log2_true: dict[str, float] = field(default_factory=dict)


def _make_species_pool(species: str, cfg: Config, rng: np.random.Generator) -> SpeciesPool:
    m = cfg["mock"]
    code = species_code(species)
    L = m["gene_length"]
    core = [(f"{code}_core_{j:03d}", random_dna(rng, L)) for j in range(m["n_core_genes"])]
    # Two named housekeeping genes carry the fluoroquinolone resistance mutations.
    core[10] = ("gyrA", random_dna(rng, L))
    core[55] = ("parC", random_dna(rng, L))
    accessory = [(f"{code}_acc_{j:03d}", random_dna(rng, L)) for j in range(m["n_accessory_genes"])]
    mutant_base = {}
    for fam, seq in (core[10], core[55]):
        mutant_base[fam] = [b for b in "ACGT" if b != seq[DETERMINANT_POS]][rng.integers(0, 3)]
    baseline = {}
    for ab in cfg["antibiotics"]:
        bp = cfg.breakpoint_log2(species, ab)
        baseline[ab] = bp - rng.integers(2, 5)          # wild type sits 2-4 dilutions below breakpoint
        if ab in m.get("intrinsic", {}).get(species, []):
            baseline[ab] += INTRINSIC_EFFECT
    return SpeciesPool(
        name=species, code=code, core=core, accessory=accessory,
        accessory_slot=rng.integers(0, len(core), len(accessory)),
        spacers=[random_dna(rng, 40) for _ in core],
        mutant_base=mutant_base, baseline=baseline,
    )


def _make_lineages(pool: SpeciesPool, cfg: Config, rng: np.random.Generator) -> list[Lineage]:
    m = cfg["mock"]
    lineages = []
    for _ in range(m["n_lineages"]):
        gyr = rng.random() < 0.35
        lineages.append(Lineage(
            accessory_profile=rng.random(len(pool.accessory)) < 0.4,
            # Each lineage tends to have or lack each plasmid -> resistance is clonal (a confounder).
            plasmid_rate={p: float(rng.choice([0.03, 0.4, 0.9])) for p in PLASMIDS},
            gyrA_fixed=gyr,
            parC_fixed=gyr and rng.random() < 0.6,
            neutral_snps=[(int(rng.integers(0, len(pool.core))), int(rng.integers(0, m["gene_length"])))
                          for _ in range(30)],
            site_bias=float(rng.beta(2, 2)),
            weight=float(rng.gamma(1.5)),
        ))
    return lineages


def _make_isolate(iid: str, pool: SpeciesPool, lineages: list[Lineage], cfg: Config,
                  shared_seqs: dict[str, str], rng: np.random.Generator) -> tuple[Isolate, str]:
    m = cfg["mock"]
    w = np.array([l.weight for l in lineages])
    li = int(rng.choice(len(lineages), p=w / w.sum()))
    lin = lineages[li]
    iso = Isolate(iid, pool.name, li)

    # --- chromosome: core genes (+ mutations) with accessory genes inserted at fixed slots
    core = [list(x) for x in pool.core]
    for gi, pos in lin.neutral_snps:                       # lineage SNPs: same base in every member
        if pos != DETERMINANT_POS:
            wt = core[gi][1][pos]
            alt = next(b for b in "ACGT"[(gi + pos) % 4:] + "ACGT" if b != wt)
            core[gi][1] = point_mutate(core[gi][1], pos, new_base=alt)
    for _ in range(5):                                     # private (isolate-specific) SNPs
        gi, pos = int(rng.integers(0, len(core))), int(rng.integers(0, m["gene_length"]))
        if pos != DETERMINANT_POS:
            core[gi][1] = point_mutate(core[gi][1], pos, rng)
    for fam, gi, fixed in (("gyrA", 10, lin.gyrA_fixed), ("parC", 55, lin.parC_fixed)):
        if rng.random() < (0.95 if fixed else 0.03):
            core[gi][1] = point_mutate(core[gi][1], DETERMINANT_POS, new_base=pool.mutant_base[fam])
            iso.carried.add(fam)

    acc_present = lin.accessory_profile ^ (rng.random(len(pool.accessory)) < 0.03)
    chromosome: list[tuple[str, str]] = []
    for ci, (fam, seq) in enumerate(core):
        chromosome.append((fam, seq + pool.spacers[ci]))
        for aj in np.flatnonzero((pool.accessory_slot == ci) & acc_present):
            chromosome.append(pool.accessory[aj])

    # --- split chromosome into contigs (as a short-read assembler would)
    n_contigs = int(rng.integers(2, 6))
    cuts = sorted(rng.choice(np.arange(1, len(chromosome)), n_contigs - 1, replace=False))
    for k, (a, b) in enumerate(zip([0, *cuts], [*cuts, len(chromosome)])):
        iso.contigs[f"{iid}_contig{k + 1}"] = chromosome[a:b]

    # --- plasmids (shared sequences across species = horizontal gene transfer)
    for pname, spec in PLASMIDS.items():
        if rng.random() < lin.plasmid_rate[pname]:
            genes = [(g, shared_seqs[g]) for g in spec["backbone"]]
            for g in spec["cargo"]:
                if rng.random() < 0.85:
                    genes.append((g, shared_seqs[g]))
                    iso.carried.add(g)
            iso.contigs[f"{iid}_{pname}"] = genes

    # --- phenotype: wild-type baseline + additive determinant effects + lab noise
    for ab in cfg["antibiotics"]:
        effect = sum(e for a, g, _, e in DETERMINANTS if a == ab and g in iso.carried)
        iso.log2_true[ab] = pool.baseline[ab] + effect + rng.normal(0, m["mic_noise_sd"])

    site = cfg["mock"]["sites"][0] if rng.random() < lin.site_bias else cfg["mock"]["sites"][1]
    return iso, site


def _format_mic(log2_value: float, lo: int, hi: int) -> str:
    v = int(np.round(log2_value))
    fmt = lambda x: f"{2.0 ** x:g}"
    if v <= lo:
        return f"<={fmt(lo)}"
    if v > hi:
        return f">{fmt(hi)}"
    return fmt(v)


def generate(cfg: Config, out: Path | None = None) -> Path:
    rng = np.random.default_rng(cfg["seed"])
    out = out or cfg.path("data_root")
    for sub in ("assemblies", "pangenome", "annotation", "reference", "truth"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    m = cfg["mock"]
    lo, hi = cfg["mic_panel"]["min_log2"], cfg["mic_panel"]["max_log2"]

    shared = sorted({g for p in PLASMIDS.values() for g in p["backbone"] + p["cargo"]})
    shared_seqs = {g: random_dna(rng, m["gene_length"]) for g in shared}

    meta_rows, mic_rows, call_rows, amr_rows = [], [], [], []
    presence: dict[str, set[str]] = {}
    n = 0
    for species in cfg["species"]:
        pool = _make_species_pool(species, cfg, rng)
        lineages = _make_lineages(pool, cfg, rng)
        for _ in range(m["isolates_per_species"]):
            n += 1
            iid = f"AMRB{n:05d}"
            iso, site = _make_isolate(iid, pool, lineages, cfg, shared_seqs, rng)

            # assemblies/<id>.fasta  +  gene_calls.tsv
            with open(out / "assemblies" / f"{iid}.fasta", "w") as fh:
                for contig, genes in iso.contigs.items():
                    seq, pos = [], 0
                    for fam, s in genes:
                        gene_len = m["gene_length"]
                        call_rows.append((iid, contig, pos + 1, pos + gene_len, "+", fam))
                        seq.append(s)
                        pos += len(s)
                    full = "".join(seq)
                    fh.write(f">{contig} len={len(full)}\n")
                    for i in range(0, len(full), 80):
                        fh.write(full[i:i + 80] + "\n")
            presence[iid] = {fam for genes in iso.contigs.values() for fam, _ in genes}

            for det in sorted(iso.carried):
                subtype = "POINT" if det in ("gyrA", "parC") else "AMR"
                amr_rows.append((iid, det, subtype))

            qc_fail = rng.random() < m["qc_fail_rate"]
            date = pd.Timestamp("2026-06-01") + pd.Timedelta(days=int(rng.integers(0, 365)))
            meta_rows.append({
                "isolate_id": iid, "species": species, "collection_site": site,
                "collection_date": date.date().isoformat(),
                "sequencing_batch": f"B{min(4, 1 + (date.month - 6) % 12 // 3)}",
                "coverage": round(float(rng.uniform(12, 25) if qc_fail else rng.uniform(35, 120)), 1),
                "n50": int(rng.uniform(8e3, 3e4) if qc_fail else rng.uniform(6e4, 3e5)),
                "contamination_pct": round(float(rng.uniform(0, 1.5)), 2),
                "qc_pass": not qc_fail,
                "_lineage": f"{pool.code}_L{iso.lineage}",
            })
            for ab, v in iso.log2_true.items():
                if rng.random() >= m["missing_rate"]:
                    mic_rows.append((iid, ab, _format_mic(v, lo, hi)))

    meta = pd.DataFrame(meta_rows)
    meta.drop(columns="_lineage").to_csv(out / "metadata.tsv", sep="\t", index=False)
    pd.DataFrame(mic_rows, columns=["isolate_id", "antibiotic", "mic"]).to_csv(out / "mic.tsv", sep="\t", index=False)
    pd.DataFrame(call_rows, columns=["isolate_id", "contig", "start", "end", "strand", "gene_family"]) \
        .to_csv(out / "pangenome" / "gene_calls.tsv", sep="\t", index=False)
    pd.DataFrame(amr_rows, columns=["isolate_id", "gene_family", "element_subtype"]) \
        .to_csv(out / "annotation" / "amrfinder_calls.tsv", sep="\t", index=False)

    families = sorted(set().union(*presence.values()))
    ids = list(presence)
    rtab = pd.DataFrame(0, index=families, columns=ids, dtype=np.uint8)
    for iid, fams in presence.items():
        rtab.loc[list(fams), iid] = 1
    rtab.index.name = "Gene"
    rtab.to_csv(out / "pangenome" / "gene_presence_absence.Rtab", sep="\t")

    pd.DataFrame([d[:3] for d in DETERMINANTS], columns=["antibiotic", "gene_family", "mechanism"]) \
        .to_csv(out / "reference" / "known_determinants.tsv", sep="\t", index=False)

    # truth/ exists ONLY for mock data: lets us check the benchmark against what we planted.
    meta[["isolate_id", "_lineage"]].rename(columns={"_lineage": "lineage"}) \
        .to_csv(out / "truth" / "lineages.tsv", sep="\t", index=False)
    decoys = sorted({g for p in PLASMIDS.values() for g in p["backbone"]})
    with open(out / "truth" / "ground_truth.json", "w") as fh:
        json.dump({
            "determinants": [dict(zip(["antibiotic", "gene_family", "mechanism", "effect_log2"], d)) for d in DETERMINANTS],
            "decoys_plasmid_backbone": decoys,
            "intrinsic": m.get("intrinsic", {}),
        }, fh, indent=2)
    return out
