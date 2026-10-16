"""The DATA CONTRACT: what the sequencing/QC pipeline must deliver, and how we read it.

This module is the seam between "data generation" and "benchmark". The mock
generator writes exactly these files; the real pipeline must write them too
(or a small adapter must convert its outputs into them).

    <data_root>/
      metadata.tsv                        one row per isolate (UCL metadata + QC report)
      mic.tsv                             long table: isolate x antibiotic -> MIC as reported by the lab
      assemblies/<isolate_id>.fasta       assembled contigs (Unicycler)
      pangenome/gene_presence_absence.Rtab  gene family x isolate 0/1 matrix (Panaroo format)
      pangenome/gene_calls.tsv            gene coordinates per isolate (Bakta GFF + Panaroo family IDs)
      annotation/amrfinder_calls.tsv      known AMR determinants found per isolate (AMRFinderPlus)
      reference/known_determinants.tsv    which gene/mutation is known to affect which antibiotic (CARD/AMRFinder DB)
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

METADATA_COLUMNS = [
    "isolate_id", "species", "collection_site", "collection_date",
    "sequencing_batch", "coverage", "n50", "contamination_pct", "qc_pass",
]
MIC_COLUMNS = ["isolate_id", "antibiotic", "mic"]
GENE_CALL_COLUMNS = ["isolate_id", "contig", "start", "end", "strand", "gene_family"]
AMRFINDER_COLUMNS = ["isolate_id", "gene_family", "element_subtype"]
KNOWN_DET_COLUMNS = ["antibiotic", "gene_family", "mechanism"]


class ContractError(ValueError):
    """Raised when input files do not match the data contract."""


def _require(df: pd.DataFrame, cols: list[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ContractError(f"{name} is missing required columns: {missing}")


# ----------------------------------------------------------------------------- readers

def load_metadata(root: Path) -> pd.DataFrame:
    df = pd.read_csv(root / "metadata.tsv", sep="\t")
    _require(df, METADATA_COLUMNS, "metadata.tsv")
    if df.isolate_id.duplicated().any():
        raise ContractError("metadata.tsv has duplicate isolate_id values")
    df["qc_pass"] = df["qc_pass"].astype(str).str.lower().isin(["true", "1", "yes"])
    return df


def parse_mic(value) -> tuple[float, str]:
    """Parse a lab MIC string into (log2 value, censor flag).

    '0.25'     -> (-2.0, '=')     exact
    '<=0.0625' -> (-4.0, '<=')    at or below the lowest tested concentration
    '>64'      -> ( 6.0, '>')     grew at the highest tested concentration
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return (np.nan, "")
    s = str(value).strip().replace(" ", "").replace("≤", "<=").replace("≥", ">=")
    prefix = next((p for p in ("<=", ">=", "<", ">") if s.startswith(p)), "")
    log2 = math.log2(float(s[len(prefix):]))
    if abs(log2 - round(log2)) > 1e-6:
        raise ContractError(f"MIC {value!r} is not a two-fold dilution step")
    v = float(round(log2))
    # Normalise everything to two censor types on the integer dilution scale:
    #   '<=v' : true MIC <= 2^v        '>v' : true MIC > 2^v (i.e. >= 2^(v+1))
    if prefix == "<":
        return (v - 1, "<=")
    if prefix == ">=":
        return (v - 1, ">")
    return (v, prefix or "=")


def load_mic(root: Path) -> pd.DataFrame:
    df = pd.read_csv(root / "mic.tsv", sep="\t", dtype={"mic": str})
    _require(df, MIC_COLUMNS, "mic.tsv")
    parsed = df["mic"].map(parse_mic)
    df["log2_mic"] = [p[0] for p in parsed]
    df["censor"] = [p[1] for p in parsed]
    return df.dropna(subset=["log2_mic"])


def load_gene_presence(root: Path) -> pd.DataFrame:
    """Panaroo/Roary Rtab -> DataFrame (rows = isolates, columns = gene families, values 0/1)."""
    rtab = pd.read_csv(root / "pangenome" / "gene_presence_absence.Rtab", sep="\t", index_col=0)
    return rtab.T.astype(np.uint8)


def load_gene_calls(root: Path) -> pd.DataFrame:
    df = pd.read_csv(root / "pangenome" / "gene_calls.tsv", sep="\t")
    _require(df, GENE_CALL_COLUMNS, "gene_calls.tsv")
    return df


def load_amrfinder(root: Path) -> pd.DataFrame:
    df = pd.read_csv(root / "annotation" / "amrfinder_calls.tsv", sep="\t")
    _require(df, AMRFINDER_COLUMNS, "amrfinder_calls.tsv")
    return df


def load_known_determinants(root: Path) -> pd.DataFrame:
    df = pd.read_csv(root / "reference" / "known_determinants.tsv", sep="\t")
    _require(df, KNOWN_DET_COLUMNS, "known_determinants.tsv")
    return df


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line.upper())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def assembly_path(root: Path, isolate_id: str) -> Path:
    return root / "assemblies" / f"{isolate_id}.fasta"


# ----------------------------------------------------------------------------- validation

def validate(root: Path, cfg) -> list[str]:
    """Check the inputs against the contract. Returns a list of human-readable findings.

    Raises ContractError for problems that make the benchmark impossible to build;
    returns warnings for things that are worth knowing but not fatal.
    """
    notes: list[str] = []
    meta = load_metadata(root)
    mic = load_mic(root)
    notes.append(f"metadata: {len(meta)} isolates, {meta.qc_pass.sum()} pass QC")

    unknown_species = set(meta.species) - set(cfg["species"])
    if unknown_species:
        notes.append(f"WARNING species not in config (will be ignored): {sorted(unknown_species)}")

    unknown_ab = set(mic.antibiotic) - set(cfg["antibiotics"])
    if unknown_ab:
        notes.append(f"WARNING antibiotics not in config (will be ignored): {sorted(unknown_ab)}")

    orphan = set(mic.isolate_id) - set(meta.isolate_id)
    if orphan:
        raise ContractError(f"{len(orphan)} isolates in mic.tsv have no metadata row")

    missing_asm = [i for i in meta.isolate_id if not assembly_path(root, i).exists()]
    if missing_asm:
        raise ContractError(f"{len(missing_asm)} isolates have no assembly, e.g. {missing_asm[:3]}")

    pa = load_gene_presence(root)
    missing_pa = set(meta.isolate_id[meta.qc_pass]) - set(pa.index)
    if missing_pa:
        raise ContractError(f"{len(missing_pa)} QC-passed isolates missing from the pan-genome matrix")

    # Coverage of the isolate x antibiotic grid (missing tests are normal; report them).
    grid = len(meta) * len(cfg["antibiotics"])
    notes.append(f"mic: {len(mic)} measurements = {len(mic) / grid:.1%} of the isolate x antibiotic grid")
    cens = mic.censor.value_counts().to_dict()
    notes.append(f"mic censoring: {cens}")
    per_ab = mic.groupby("antibiotic").size().reindex(cfg["antibiotics"]).fillna(0).astype(int)
    notes.append("mic per antibiotic: " + ", ".join(f"{k}={v}" for k, v in per_ab.items()))
    return notes
