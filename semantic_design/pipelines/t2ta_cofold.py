"""
Run ESMFold multimer + pDockQ scoring on toxin-antitoxin protein pairs.

Usage: python pipelines/t2ta_cofold.py --config path/to/config.yaml
"""

import argparse
import hashlib
import os
import sys
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Tuple

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

import numpy as np
import pandas as pd
import torch
import yaml
from tqdm import tqdm

import esm


def sanitize_identifier(value: str) -> str:
    """Return a filesystem-friendly identifier derived from arbitrary text.

    Args:
        value: Arbitrary ID string.

    Returns:
        Identifier containing alphanumerics/-/_ suitable for filenames.
    """
    safe = "".join(ch for ch in value if ch.isalnum() or ch in ("-", "_"))
    return safe or "pair"


def shorten_identifier(value: str, max_length: int = 120) -> str:
    """Ensure identifiers remain within filesystem limits using hashing."""
    sanitized = sanitize_identifier(value)
    if len(sanitized) <= max_length:
        return sanitized
    digest = hashlib.md5(sanitized.encode("utf-8")).hexdigest()[:10]
    keep = max_length - len(digest) - 1
    return f"{sanitized[:keep]}_{digest}"


@dataclass
class CofoldConfig:
    """Configuration for toxin–antitoxin co-folding + pDockQ scoring."""

    pairs_csv: Path  # CSV produced by identify_unique_pairs
    output_dir: Path  # Directory for co-fold inputs, PDBs, and score summaries

    # Column names used in the pairs CSV
    root_id_col: str = "Root_ID"
    sequence1_col: str = "Amino Acid Sequence 1"
    sequence2_col: str = "Amino Acid Sequence 2"
    sequence1_id_col: str = "Evo Sequence ID 1"
    sequence2_id_col: str = "Evo Sequence ID 2"

    run_esmfold: bool = True
    run_pdockq: bool = True
    pdockq_threshold: float = 0.23

    cofold_input_csv: Path = field(init=False)
    esmfold_output_dir: Path = field(init=False)
    pdockq_scores_csv: Path = field(init=False)
    pdockq_summary_csv: Path = field(init=False)
    pdockq_filtered_csv: Path = field(init=False)
    pdockq_filtered_fasta: Path = field(init=False)

    def __post_init__(self) -> None:
        """Resolve key directories and set every downstream output path."""
        self.pairs_csv = Path(self.pairs_csv)
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.cofold_input_csv = self.output_dir / "cofold_input.csv"
        self.esmfold_output_dir = self.output_dir / "esmfold_structures"
        self.esmfold_output_dir.mkdir(parents=True, exist_ok=True)
        self.pdockq_scores_csv = self.output_dir / "pdockq_scores.csv"
        self.pdockq_summary_csv = self.output_dir / "pdockq_summary.csv"
        self.pdockq_filtered_csv = self.output_dir / "pdockq_high_confidence.csv"
        self.pdockq_filtered_fasta = self.output_dir / "pdockq_high_confidence.fasta"


def load_config(config_path: Path) -> CofoldConfig:
    """Load YAML configuration file."""
    with open(config_path, "r") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Configuration must be a mapping: {config_path}")
    return CofoldConfig(**data)


def prepare_cofold_inputs(config: CofoldConfig) -> pd.DataFrame:
    """Convert toxin-antitoxin pair CSV into esmfold_multimer input format.

    Args:
        config: Cofold configuration describing column names and output paths.

    Returns:
        DataFrame with `sequence1`, `sequence2`, and `id_pair` columns used for folding.
    """
    df = pd.read_csv(config.pairs_csv)
    required_columns = [
        config.root_id_col,
        config.sequence1_col,
        config.sequence2_col,
        config.sequence1_id_col,
        config.sequence2_id_col,
    ]
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {config.pairs_csv}: {missing}")

    def build_id(row: pd.Series) -> str:
        seq1_id = str(row[config.sequence1_id_col])
        seq2_id = str(row[config.sequence2_id_col])
        raw_id = f"{seq1_id}_{seq2_id}"
        return shorten_identifier(raw_id)

    cofold_df = pd.DataFrame(
        {
            "sequence1": df[config.sequence1_col].astype(str),
            "sequence2": df[config.sequence2_col].astype(str),
            "id_pair": [build_id(row) for _, row in df.iterrows()],
        }
    )
    cofold_df = cofold_df.drop_duplicates(subset="id_pair").reset_index(drop=True)
    cofold_df.to_csv(config.cofold_input_csv, index=False)
    return cofold_df


def run_esmfold(config: CofoldConfig, cofold_df: pd.DataFrame) -> None:
    """Fold toxin/antitoxin pairs using the bundled ESMFold helper."""

    MIN_SEQUENCE_LENGTH = 100
    MAX_SEQUENCE_LENGTH = 1024

    def fold_multimer(seq1: str, seq2: str, pair_id: str, model: "ESMFold") -> None:
        total_len = len(seq1) + len(seq2)
        if total_len < MIN_SEQUENCE_LENGTH or total_len > MAX_SEQUENCE_LENGTH:
            return

        pair_dir = config.esmfold_output_dir / pair_id
        pdb_path = pair_dir / f"{pair_id}.pdb"
        if pdb_path.exists():
            return

        multimer_seq = f"{seq1}:{seq2}"
        with torch.no_grad():
            pdb = model.infer_pdb(multimer_seq)

        pair_dir.mkdir(parents=True, exist_ok=True)
        with open(pdb_path, "w") as handle:
            handle.write(pdb)

    model = esm.pretrained.esmfold_v1()
    model = model.eval().to("cuda:0")
    model.set_chunk_size(128)
    model.half()

    for seq1, seq2, pair_id in tqdm(
        list(zip(cofold_df["sequence1"], cofold_df["sequence2"], cofold_df["id_pair"])),
        desc="Folding pairs",
    ):
        fold_multimer(seq1, seq2, pair_id, model)


def extract_pdockq_scores(config: CofoldConfig) -> None:
    """Compute pDockQ metrics for every generated PDB and write a CSV."""

    def parse_atm_record(line: str) -> Dict[str, Any]:
        record: Dict[str, Any] = defaultdict()
        record["atm_name"] = line[12:16].strip()
        record["res_name"] = line[17:20].strip()
        record["chain"] = line[21]
        record["res_no"] = int(line[22:26])
        record["coords"] = (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )
        record["B"] = float(line[60:66])
        return record

    def pdb_to_coords(pdb: str) -> Tuple[Dict[str, List[List[float]]], np.ndarray]:
        chain_coords: Dict[str, List[List[float]]] = defaultdict(list)
        plddt_dict: "OrderedDict[str, List[float]]" = OrderedDict()
        for line in pdb.splitlines():
            if not line.startswith("ATOM"):
                continue
            rec = parse_atm_record(line)
            if rec["atm_name"] == "CB" or (rec["atm_name"] == "CA" and rec["res_name"] == "GLY"):
                chain_coords[rec["chain"]].append(list(rec["coords"]))
                res_id = f"{rec['chain']}{rec['res_no']}"
                plddt_dict.setdefault(res_id, []).append(rec["B"])
        plddt = np.array([np.mean(vals) for vals in plddt_dict.values()])
        return chain_coords, plddt

    def calc_pdockq(chain_coords: Dict[str, List[List[float]]], plddt: np.ndarray) -> Tuple[float, float, int, float]:
        chains = list(chain_coords.keys())
        if len(chains) < 2 or plddt.size == 0:
            return 0.0, 0.0, 0, 0.0
        coords1 = np.array(chain_coords[chains[0]])
        coords2 = np.array(chain_coords[chains[1]])
        mat = np.append(coords1, coords2, axis=0)
        diffs = mat[:, np.newaxis, :] - mat[np.newaxis, :, :]
        dists = np.sqrt(np.sum(diffs**2, axis=-1))
        l1 = len(coords1)
        contact_dists = dists[:l1, l1:]
        contacts = np.argwhere(contact_dists <= 8)
        if contacts.size == 0:
            return 0.0, 0.0, 0, float(plddt.mean())
        avg_if_plddt = float(
            np.average(
                np.concatenate(
                    [plddt[np.unique(contacts[:, 0])], plddt[np.unique(contacts[:, 1])]]
                )
            )
        )
        n_if_contacts = int(contacts.shape[0])
        x = avg_if_plddt * np.log10(n_if_contacts + 1)
        pdockq = float(0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018)
        return pdockq, avg_if_plddt, n_if_contacts, float(plddt.mean())

    records: List[Dict[str, Any]] = []
    for folder_name, _, filenames in os.walk(config.esmfold_output_dir):
        for filename in filenames:
            if not filename.endswith(".pdb"):
                continue
            pdb_path = Path(folder_name) / filename
            with open(pdb_path, "r") as handle:
                pdb_contents = handle.read()
            chain_coords, plddt = pdb_to_coords(pdb_contents)
            pdockq, if_plddt, if_contacts, avg_plddt = calc_pdockq(chain_coords, plddt)
            records.append(
                {
                    "PDB_File": str(pdb_path),
                    "pDockQ": pdockq,
                    "if_pLDDT": if_plddt,
                    "if_contacts": if_contacts,
                    "avg_pLDDT": avg_plddt,
                }
            )

    pd.DataFrame(records).to_csv(config.pdockq_scores_csv, index=False)


def summarize_pdockq(
    config: CofoldConfig, cofold_df: pd.DataFrame, pdockq_df: pd.DataFrame
) -> None:
    """Merge pDockQ scores with sequence metadata and emit CSV/FASTA reports.

    Args:
        config: Cofold configuration with score thresholds and paths.
        cofold_df: Precomputed toxin/antitoxin sequence table.
        pdockq_df: Raw pDockQ metrics DataFrame.

    Returns:
        None.
    """
    if pdockq_df.empty:
        print("No pDockQ scores were extracted. Skipping summary generation.")
        return

    pdockq_df["id_pair"] = pdockq_df["PDB_File"].apply(lambda p: Path(p).stem)
    for col in ["pDockQ", "if_pLDDT", "if_contacts", "avg_pLDDT"]:
        pdockq_df[col] = pd.to_numeric(pdockq_df[col], errors="coerce")

    merged = pdockq_df.merge(cofold_df, on="id_pair", how="left")
    merged = merged.sort_values(by="pDockQ", ascending=False)
    merged.to_csv(config.pdockq_summary_csv, index=False)

    high_conf = merged[merged["pDockQ"] >= config.pdockq_threshold]
    high_conf.to_csv(config.pdockq_filtered_csv, index=False)

    if high_conf.empty:
        print(
            f"No complexes exceeded pDockQ >= {config.pdockq_threshold}. "
            "Skipping FASTA export."
        )
        return

    with open(config.pdockq_filtered_fasta, "w") as fasta_handle:
        for _, row in high_conf.iterrows():
            fasta_handle.write(f">{row['id_pair']}_tox\n{row['sequence1']}\n")
            fasta_handle.write(f">{row['id_pair']}_antitox\n{row['sequence2']}\n")


def run_pipeline(config: CofoldConfig) -> None:
    """Main entry point for co-fold execution, pDockQ scoring, and reporting.

    Args:
        config: Parsed `CofoldConfig` describing all IO locations.

    Returns:
        None. Outputs reside in `config.output_dir`.
    """
    cofold_df = prepare_cofold_inputs(config)

    if config.run_esmfold:
        print("Running ESMFold multimer on toxin-antitoxin pairs...", flush=True)
        run_esmfold(config, cofold_df)
    else:
        print("Skipping ESMFold run (run_esmfold=False).", flush=True)

    if config.run_pdockq:
        print("Extracting pDockQ metrics...", flush=True)
        extract_pdockq_scores(config)
        pdockq_df = pd.read_csv(config.pdockq_scores_csv)
        summarize_pdockq(config, cofold_df, pdockq_df)
    else:
        print("Skipping pDockQ extraction (run_pdockq=False).", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run ESMFold multimer + pDockQ scoring on toxin-antitoxin pairs."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML config describing pair CSV and ESMFold repo location.",
    )
    args = parser.parse_args()
    run_pipeline(load_config(Path(args.config)))
