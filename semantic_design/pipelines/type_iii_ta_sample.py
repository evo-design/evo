"""
Type III toxin-antitoxin sampling pipeline using Evo.

Usage:
    python pipelines/type_iii_ta_sample.py --config path/to/config.yaml

This pipeline mirrors the sampling/filtering flow of the other semantic design
pipelines and adds Type III-specific analyses (tandem repeat finder, RNA folding,
and hairpin-based pairing heuristics).
"""

import argparse
import ast
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import ViennaRNA
import yaml
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

DEFAULT_RNA_STRUCTURE_FILTER_SCRIPT = Path(
    "/large_storage/hielab/userspace/adititm/semantic_mining_data/toxin_antitoxin/type_iii/rebuttal/find_matching_rna_strucs.py"
)
DEFAULT_RNA_SEQUENCE_FILTER_SCRIPT = Path(
    "/large_storage/hielab/userspace/adititm/semantic_mining_data/toxin_antitoxin/type_iii/rebuttal/find_matching_rna_seqs.py"
)
DEFAULT_PFAM_REFERENCE_CSV = Path(
    "/large_storage/hielab/userspace/adititm/semantic_mining_data/toxin_antitoxin/type_iii/wt/type_III_T_wt_hmm_hits.csv"
)

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from semantic_design import (
    read_prompts,
    model_load,
    sample_model,
    get_rc,
    make_fasta,
    run_prodigal,
    filter_protein_fasta,
    fold_proteins,
    filter_proteins_by_threshold,
)


@dataclass
class Config:
    """Configuration for Type III toxin-antitoxin generation and downstream analysis."""

    input_prompts: Path
    output_dir: Path
    segmasker_path: Path
    trf_path: Path  # Path to tandem repeat finder binary
    rna_structures_reference_csv: Path  # CSV with experimental RNA folds for comparison

    model_name: str
    n_tokens: int
    temperature: float
    top_k: int
    batched: bool
    batch_size: int
    n_sample_per_prompt: int

    rc_truth: bool
    return_both: bool
    filter_min_length: int = 50
    filter_max_length: int = 400
    filter_partial_bool: bool = False
    segmasker_threshold: float = 0.2
    run_esm_fold: bool = True
    plddt_threshold: float = 0.3
    ptm_threshold: float = 0.0

    write_trf_to_csv: bool = True
    rna_structure_filter_reference_csv: Optional[Path] = None
    rna_structure_filter_script: Path = DEFAULT_RNA_STRUCTURE_FILTER_SCRIPT
    rna_structure_filter_structure_type: str = "both"
    rna_structure_filter_min_similarity: float = 0.7
    rna_structure_filter_pre_filter_threshold: float = 0.7
    rna_structure_filter_batch_size: int = 100
    rna_structure_filter_max_results: Optional[int] = None
    rna_structure_filter_cpus: Optional[int] = None
    rna_sequence_filter_reference_csv: Optional[Path] = None
    rna_sequence_filter_script: Path = DEFAULT_RNA_SEQUENCE_FILTER_SCRIPT
    rna_sequence_filter_min_identity: float = 70.0
    rna_sequence_filter_processes: Optional[int] = None
    hmmscan_binary: str = "hmmscan"
    hmmscan_pfam_db_path: Optional[Path] = None
    hmmscan_cpu: int = 4
    pfam_reference_hits_csv: Path = DEFAULT_PFAM_REFERENCE_CSV
    pfam_evalue_threshold: float = 0.05
    rna_require_hairpin: bool = True
    rna_minimum_mfe: float = -3.0
    rna_require_all_bases: bool = True
    cmscan_binary: str = "cmscan"
    cmscan_model_paths: Optional[List[Path]] = None
    cmscan_evalue_threshold: float = 0.05
    cmscan_allowed_families: Optional[List[str]] = None
    cmscan_allowed_families_csv: Optional[Path] = None
    cmscan_allowed_families_column: str = "Query Name"

    evo_gen_seqs_file_save_location: Path = field(init=False)
    all_seqs_fasta: Path = field(init=False)
    proteins_file: Path = field(init=False)
    orfs_file: Path = field(init=False)
    filtered_proteins_file: Path = field(init=False)
    output_folds_file: Path = field(init=False)
    output_filtered_folds: Path = field(init=False)
    output_trf_csv: Path = field(init=False)
    rna_fold_csv: Path = field(init=False)
    ta_pairs_csv: Path = field(init=False)
    rna_candidates_csv: Path = field(init=False)
    rna_structure_matches_csv: Path = field(init=False)
    rna_sequence_matches_csv: Path = field(init=False)
    rna_candidates_fasta: Path = field(init=False)
    hmmscan_domtblout: Path = field(init=False)
    hmmscan_hits_csv: Path = field(init=False)
    cmscan_tblout_dir: Path = field(init=False)
    cmscan_hits_csv: Path = field(init=False)
    final_candidates_csv: Path = field(init=False)

    def __post_init__(self) -> None:
        self.input_prompts = Path(self.input_prompts)
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.segmasker_path = Path(self.segmasker_path)
        self.trf_path = Path(self.trf_path)
        self.rna_structures_reference_csv = Path(self.rna_structures_reference_csv)
        if self.rna_structure_filter_reference_csv:
            self.rna_structure_filter_reference_csv = Path(self.rna_structure_filter_reference_csv)
        if self.rna_sequence_filter_reference_csv:
            self.rna_sequence_filter_reference_csv = Path(self.rna_sequence_filter_reference_csv)
        if self.hmmscan_pfam_db_path:
            self.hmmscan_pfam_db_path = Path(self.hmmscan_pfam_db_path)
        if self.pfam_reference_hits_csv:
            self.pfam_reference_hits_csv = Path(self.pfam_reference_hits_csv)
        if self.cmscan_model_paths:
            self.cmscan_model_paths = [Path(path) for path in self.cmscan_model_paths]
        if isinstance(self.cmscan_allowed_families, str):
            self.cmscan_allowed_families = [self.cmscan_allowed_families]
        if self.cmscan_allowed_families_csv:
            self.cmscan_allowed_families_csv = Path(self.cmscan_allowed_families_csv)

        self.evo_gen_seqs_file_save_location = self.output_dir / "generated_sequences.csv"
        self.all_seqs_fasta = self.output_dir / "all_sequences.fasta"
        self.proteins_file = self.output_dir / "proteins.fasta"
        self.orfs_file = self.output_dir / "orfs.fasta"
        self.filtered_proteins_file = self.output_dir / "filtered_proteins.fasta"
        self.output_folds_file = self.output_dir / "folds.csv"
        self.output_filtered_folds = self.output_dir / "filtered_folds.csv"
        self.output_trf_csv = self.output_dir / "tandem_repeats.csv"
        self.rna_fold_csv = self.output_dir / "rna_fold_predictions.csv"
        self.ta_pairs_csv = self.output_dir / "type_iii_pairs.csv"
        self.rna_candidates_csv = self.output_dir / "rna_candidates.csv"
        self.rna_structure_matches_csv = self.output_dir / "rna_structure_matches.csv"
        self.rna_sequence_matches_csv = self.output_dir / "rna_sequence_matches.csv"
        self.rna_candidates_fasta = self.output_dir / "rna_candidates.fasta"
        self.hmmscan_domtblout = self.output_dir / "hmmscan.domtblout"
        self.hmmscan_hits_csv = self.output_dir / "hmmscan_hits.csv"
        self.cmscan_tblout_dir = self.output_dir / "cmscan_tblout"
        self.cmscan_tblout_dir.mkdir(parents=True, exist_ok=True)
        self.cmscan_hits_csv = self.output_dir / "cmscan_hits.csv"
        self.final_candidates_csv = self.output_dir / "filtered_type_iii_candidates.csv"


def load_config(config_path: Path) -> Config:
    with open(config_path, "r") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Configuration must be a mapping: {config_path}")
    return Config(**data)


def load_generated_sequences(csv_path: Path) -> Dict[str, str]:
    df = pd.read_csv(csv_path)
    if "UUID" in df.columns and "Generated Sequence" in df.columns:
        return dict(zip(df["UUID"].astype(str), df["Generated Sequence"].astype(str)))
    # Fallback to positional columns
    return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 2].astype(str)))


def run_tandem_repeat_finder(sequence: str, root_id: str, trf_path: Path) -> pd.DataFrame:
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as temp_input_file:
        temp_input_file.write(f">sequence\n{sequence}\n")
    result = subprocess.run(
        [
            str(trf_path),
            temp_input_file.name,
            "2",
            "7",
            "7",
            "80",
            "10",
            "50",
            "500",
            "-h",
            "-ngs",
        ],
        capture_output=True,
        text=True,
    )
    repeats: List[Dict[str, Any]] = []
    for line in result.stdout.splitlines():
        if line.startswith("@"):
            continue
        data = line.strip().split()
        if len(data) < 14:
            continue
        repeats.append(
            {
                "Root ID": root_id,
                "Start": int(data[0]),
                "End": int(data[1]),
                "Period Size": float(data[2]),
                "Copy Number": float(data[3]),
                "Consensus Size": int(data[4]),
                "Percent Match": float(data[5]),
                "Percent Indels": float(data[6]),
                "Alignment Score": int(data[7]),
                "A": int(data[8]),
                "C": int(data[9]),
                "G": int(data[10]),
                "T": int(data[11]),
                "Entropy": float(data[12]),
                "Repeat Sequence": data[13],
                "Full TRF Region": sequence[int(data[0]) : int(data[1])],
            }
        )
    Path(temp_input_file.name).unlink(missing_ok=True)
    return pd.DataFrame(repeats)


def get_tandem_repeats(
    filtered_folds: pd.DataFrame,
    sequences_csv: Path,
    config: Config,
) -> pd.DataFrame:
    if filtered_folds.empty:
        return pd.DataFrame()

    seq_map = load_generated_sequences(sequences_csv)
    trf_frames = []
    for _, row in filtered_folds.iterrows():
        seq_id = str(row["Evo Sequence ID"])
        root_id = seq_id.split("_")[0]
        dna_seq = seq_map.get(root_id)
        if not dna_seq:
            continue
        trf_result = run_tandem_repeat_finder(dna_seq, root_id, config.trf_path)
        if not trf_result.empty:
            trf_frames.append(trf_result)

    if not trf_frames:
        return pd.DataFrame()

    trf_df = pd.concat(trf_frames, ignore_index=True)
    if config.write_trf_to_csv:
        trf_df.to_csv(config.output_trf_csv, index=False)
    return trf_df


def predict_rna_structure(rna_sequence: str) -> Tuple[str, float]:
    fold_compound = ViennaRNA.fold_compound(rna_sequence)
    structure, mfe = fold_compound.mfe()
    return structure, mfe


def predict_hairpins(dot_bracket: str) -> List[Tuple[int, int, int, int]]:
    hairpins: List[Tuple[int, int, int, int]] = []
    stack: List[int] = []
    for idx, char in enumerate(dot_bracket):
        if char == "(":
            stack.append(idx)
        elif char == ")" and stack:
            start = stack.pop()
            end = idx
            if end - start <= 4:
                continue
            loop_start = start + 1
            loop_end = end - 1
            if all(dot_bracket[pos] == "." for pos in range(loop_start, loop_end + 1)):
                hairpins.append((start, end, loop_start, loop_end))
    return hairpins


def fold_trfs(trf_df: pd.DataFrame, output_csv: Path) -> pd.DataFrame:
    if trf_df.empty:
        empty = pd.DataFrame(
            columns=[
                "Evo Sequence ID",
                "Description",
                "DNA Sequence",
                "RNA Sequence",
                "Secondary Structure",
                "MFE",
                "Hairpins",
            ]
        )
        empty.to_csv(output_csv, index=False)
        return empty

    records = trf_df["Full TRF Region"].astype(str).tolist()
    ids = trf_df["Root ID"].astype(str).tolist()
    fold_info: List[Dict[str, Any]] = []
    for dna_seq, seq_id in zip(records, ids):
        rna_seq = str(Seq(dna_seq).transcribe())
        structure, mfe = predict_rna_structure(rna_seq)
        hairpins = predict_hairpins(structure)
        fold_info.append(
            {
                "Evo Sequence ID": seq_id,
                "Description": seq_id,
                "DNA Sequence": dna_seq,
                "RNA Sequence": rna_seq,
                "Secondary Structure": structure,
                "MFE": mfe,
                "Hairpins": hairpins,
            }
        )
    fold_df = pd.DataFrame(fold_info)
    fold_df.to_csv(output_csv, index=False)
    return fold_df


def get_at_pairs(rna_fold_df: pd.DataFrame, filtered_folds: pd.DataFrame, output_csv: Path) -> pd.DataFrame:
    if rna_fold_df.empty or filtered_folds.empty:
        empty = pd.DataFrame()
        empty.to_csv(output_csv, index=False)
        return empty

    rna_fold_df["Has Hairpin"] = rna_fold_df["Hairpins"].apply(lambda hp: bool(hp))
    hairpin_df = rna_fold_df[rna_fold_df["Has Hairpin"]].copy()
    filtered_folds["Root ID"] = filtered_folds["Evo Sequence ID"].astype(str).str.split("_").str[0]
    merged = hairpin_df.merge(
        filtered_folds,
        left_on="Evo Sequence ID",
        right_on="Root ID",
        how="inner",
        suffixes=("", "_Protein"),
    )
    merged.to_csv(output_csv, index=False)
    return merged


def visualize_rna_structures(reference_csv: Path, rna_fold_df: pd.DataFrame) -> Tuple[List[np.ndarray], List[str]]:
    if rna_fold_df.empty:
        return [], []
    reference = pd.read_csv(reference_csv)
    combined = pd.concat([rna_fold_df, reference], ignore_index=True, sort=False)
    combined = combined[combined["RNA Sequence"].apply(lambda seq: len(str(seq)) <= 500)]

    flattened = []
    labels = []
    for _, row in combined.iterrows():
        rna_seq = str(row["RNA Sequence"])
        structure = str(row["Secondary Structure"])
        ViennaRNA.fold_compound(rna_seq)  # ensure ViennaRNA caches
        _, ensemble_energy = ViennaRNA.pf_fold(rna_seq)
        prob_matrix = np.zeros((len(rna_seq) + 1, len(rna_seq) + 1))
        for i in range(1, len(rna_seq)):
            for j in range(i + 1, len(rna_seq) + 1):
                prob_matrix[i, j] = ViennaRNA.get_pr(i, j)
        flattened.append(np.sum(prob_matrix, axis=0))

        desc = str(row.get("Description", "Evo Generated")).lower()
        if "ran" in desc:
            labels.append("Random")
        elif "exp" in desc:
            labels.append("Experimental")
        elif "hyp" in desc:
            labels.append("Hypothetical")
        else:
            labels.append("Evo Generated")

    return flattened, labels


def _hairpin_exists(hairpin_value: Any) -> bool:
    if isinstance(hairpin_value, list):
        return bool(hairpin_value)
    if isinstance(hairpin_value, str):
        value = hairpin_value.strip()
        if not value or value == "[]":
            return False
        try:
            parsed = ast.literal_eval(value)
            return bool(parsed)
        except (ValueError, SyntaxError):
            return False
    return bool(hairpin_value)


def _has_all_bases(dna_sequence: Any) -> bool:
    seq = str(dna_sequence or "").upper()
    return all(base in seq for base in "ACGT")


def filter_folded_trfs(
    trf_df: pd.DataFrame,
    fold_trf_df: pd.DataFrame,
    config: Config,
) -> Tuple[pd.DataFrame, Set[str]]:
    if fold_trf_df.empty:
        return fold_trf_df, set()

    mask = pd.Series(True, index=fold_trf_df.index)
    if config.rna_require_hairpin:
        mask &= fold_trf_df["Hairpins"].apply(_hairpin_exists)
    if config.rna_minimum_mfe is not None:
        mask &= fold_trf_df["MFE"].apply(lambda val: pd.notna(val) and float(val) <= config.rna_minimum_mfe)
    if config.rna_require_all_bases:
        mask &= fold_trf_df["DNA Sequence"].apply(_has_all_bases)

    filtered = fold_trf_df[mask].copy()
    filtered.to_csv(config.rna_fold_csv, index=False)
    passing_ids: Set[str] = set(filtered["Evo Sequence ID"].astype(str))
    return filtered, passing_ids


def write_rna_candidates_fasta(candidate_table: pd.DataFrame, fasta_path: Path) -> None:
    if candidate_table.empty:
        if fasta_path.exists():
            fasta_path.unlink()
        return

    records: List[SeqRecord] = []
    for _, row in candidate_table.iterrows():
        rna_seq = str(row.get("Full RNA", "")).replace("U", "U")
        if not rna_seq:
            continue
        seq_id = str(row.get("Sequence_ID", row.get("Root ID", "candidate")))
        description = str(row.get("Root ID", seq_id))
        records.append(SeqRecord(Seq(rna_seq), id=seq_id, description=description))

    if records:
        SeqIO.write(records, fasta_path, "fasta")
    else:
        if fasta_path.exists():
            fasta_path.unlink()


def prepare_rna_candidate_table(trf_df: pd.DataFrame, fold_trf_df: pd.DataFrame, output_csv: Path) -> pd.DataFrame:
    if trf_df.empty or fold_trf_df.empty:
        empty = pd.DataFrame()
        empty.to_csv(output_csv, index=False)
        return empty

    rename_map = {
        "Period Size": "Period_size",
        "Copy Number": "Copy_number",
        "Consensus Size": "Consensus_size",
        "Percent Match": "Percent_match",
        "Percent Indels": "Percent_indels",
        "Alignment Score": "Alignment_score",
        "Repeat Sequence": "Sequence",
    }
    trf_standard = trf_df.rename(columns={k: v for k, v in rename_map.items() if k in trf_df.columns}).copy()
    if "Sequence_ID" not in trf_standard.columns:
        trf_standard["Sequence_ID"] = trf_standard.apply(
            lambda row: f"{row['Root ID']}_{int(row['Start'])}_{int(row['End'])}", axis=1
        )

    merged = trf_standard.merge(
        fold_trf_df,
        left_on="Root ID",
        right_on="Evo Sequence ID",
        how="inner",
        suffixes=("", "_RNA"),
    )

    if merged.empty:
        merged.to_csv(output_csv, index=False)
        return merged

    merged["Full RNA"] = merged["RNA Sequence"]
    merged["Full Structure"] = merged["Secondary Structure"]
    merged["Full MFE"] = merged["MFE"]
    merged["Full Hairpins"] = merged["Hairpins"]
    for consensus_col, full_col in [
        ("Consensus RNA", "Full RNA"),
        ("Consensus Structure", "Full Structure"),
        ("Consensus MFE", "Full MFE"),
        ("Consensus Hairpins", "Full Hairpins"),
    ]:
        if consensus_col not in merged.columns:
            merged[consensus_col] = merged[full_col]

    columns_order = [
        "Root ID",
        "Sequence_ID",
        "Start",
        "End",
        "Period_size",
        "Copy_number",
        "Consensus_size",
        "Percent_match",
        "Percent_indels",
        "Alignment_score",
        "A",
        "C",
        "G",
        "T",
        "Entropy",
        "Sequence",
        "Full TRF Region",
        "Consensus RNA",
        "Consensus Structure",
        "Consensus MFE",
        "Consensus Hairpins",
        "Full RNA",
        "Full Structure",
        "Full MFE",
        "Full Hairpins",
    ]
    available_columns = [col for col in columns_order if col in merged.columns]
    merged.loc[:, available_columns].to_csv(output_csv, index=False)
    return merged


def run_rna_structure_filter(candidates_csv: Path, config: Config) -> Set[str]:
    script_path = config.rna_structure_filter_script
    if not script_path or not script_path.exists():
        return set()

    target_csv = config.rna_structure_filter_reference_csv or config.rna_structures_reference_csv
    if not target_csv or not target_csv.exists():
        return set()

    if not candidates_csv.exists():
        return set()

    cmd = [
        sys.executable,
        str(script_path),
        "--query",
        str(candidates_csv),
        "--target",
        str(target_csv),
        "--output",
        str(config.rna_structure_matches_csv),
        "--structure-type",
        config.rna_structure_filter_structure_type,
        "--min-similarity",
        str(config.rna_structure_filter_min_similarity),
        "--pre-filter-threshold",
        str(config.rna_structure_filter_pre_filter_threshold),
        "--batch-size",
        str(config.rna_structure_filter_batch_size),
    ]
    if config.rna_structure_filter_max_results:
        cmd.extend(["--max-results", str(config.rna_structure_filter_max_results)])
    if config.rna_structure_filter_cpus:
        cmd.extend(["--cpus", str(config.rna_structure_filter_cpus)])

    subprocess.run(cmd, check=True)
    if not config.rna_structure_matches_csv.exists():
        return set()

    results = pd.read_csv(config.rna_structure_matches_csv)
    if "Query_ID" not in results.columns:
        return set()
    return set(results["Query_ID"].astype(str))


def run_rna_sequence_filter(candidates_csv: Path, config: Config) -> Set[str]:
    script_path = config.rna_sequence_filter_script
    if not script_path or not script_path.exists():
        return set()

    reference_csv = config.rna_sequence_filter_reference_csv or config.rna_structures_reference_csv
    if not reference_csv or not reference_csv.exists():
        return set()

    if not candidates_csv.exists():
        return set()

    cmd = [
        sys.executable,
        str(script_path),
        "--reference_csv",
        str(reference_csv),
        "--comparison_csv",
        str(candidates_csv),
        "--output_csv",
        str(config.rna_sequence_matches_csv),
        "--min-identity",
        str(config.rna_sequence_filter_min_identity),
    ]
    if config.rna_sequence_filter_processes:
        cmd.extend(["--processes", str(config.rna_sequence_filter_processes)])

    subprocess.run(cmd, check=True)
    if not config.rna_sequence_matches_csv.exists():
        return set()

    results = pd.read_csv(config.rna_sequence_matches_csv)
    if "comp_root_id" not in results.columns:
        return set()
    filtered = results[results.get("identity_percent", 0) >= config.rna_sequence_filter_min_identity]
    return set(filtered["comp_root_id"].astype(str))


def load_allowed_pfam_names(csv_path: Optional[Path]) -> Set[str]:
    if not csv_path or not csv_path.exists():
        return set()
    df = pd.read_csv(csv_path)
    if "pfam_name" not in df.columns:
        return set()
    return {str(name).strip().strip('"') for name in df["pfam_name"].dropna().unique()}


def parse_domtblout(domtbl_path: Path) -> pd.DataFrame:
    if not domtbl_path.exists():
        return pd.DataFrame()
    hits: List[Dict[str, Any]] = []
    with open(domtbl_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 23:
                continue
            hit = {
                "pfam_id": fields[0],
                "pfam_accession": fields[1],
                "sequence_id": fields[3],
                "e_value": float(fields[11]) if fields[11] != "-" else None,
                "score": float(fields[13]) if fields[13] != "-" else None,
                "bias": float(fields[14]) if fields[14] != "-" else None,
                "hmm_from": int(fields[15]) if fields[15] != "-" else None,
                "hmm_to": int(fields[16]) if fields[16] != "-" else None,
                "ali_from": int(fields[17]) if fields[17] != "-" else None,
                "ali_to": int(fields[18]) if fields[18] != "-" else None,
                "pfam_name": " ".join(fields[22:]) if len(fields) > 22 else "",
            }
            hits.append(hit)
    return pd.DataFrame(hits)


def run_hmmscan_filter(config: Config) -> pd.DataFrame:
    if not config.hmmscan_pfam_db_path or not Path(config.hmmscan_pfam_db_path).exists():
        return pd.DataFrame()
    if not Path(config.filtered_proteins_file).exists():
        return pd.DataFrame()

    cmd = [
        config.hmmscan_binary,
        "--domtblout",
        str(config.hmmscan_domtblout),
        "--cpu",
        str(config.hmmscan_cpu),
        str(config.hmmscan_pfam_db_path),
        str(config.filtered_proteins_file),
    ]
    subprocess.run(cmd, check=True)
    hits_df = parse_domtblout(config.hmmscan_domtblout)
    hits_df.to_csv(config.hmmscan_hits_csv, index=False)

    allowed_names = load_allowed_pfam_names(config.pfam_reference_hits_csv)
    if not allowed_names or hits_df.empty:
        if config.pfam_evalue_threshold is not None and not hits_df.empty:
            hits_df = hits_df[hits_df["e_value"].apply(lambda val: pd.notna(val) and float(val) <= config.pfam_evalue_threshold)]
            hits_df.to_csv(config.hmmscan_hits_csv, index=False)
        return hits_df

    def _normalize(value: Any) -> str:
        return str(value).strip().strip('"')

    normalized_allowed = {_normalize(name) for name in allowed_names}
    mask = hits_df["pfam_name"].apply(lambda name: _normalize(name) in normalized_allowed)
    filtered = hits_df[mask]
    if config.pfam_evalue_threshold is not None and not filtered.empty:
        filtered = filtered[filtered["e_value"].apply(lambda val: pd.notna(val) and float(val) <= config.pfam_evalue_threshold)]
    filtered.to_csv(config.hmmscan_hits_csv, index=False)
    return filtered


def load_allowed_families(config: Config) -> Set[str]:
    allowed: Set[str] = set()
    if config.cmscan_allowed_families:
        allowed.update({str(name).strip().lower() for name in config.cmscan_allowed_families})
    if config.cmscan_allowed_families_csv and config.cmscan_allowed_families_csv.exists():
        df = pd.read_csv(config.cmscan_allowed_families_csv)
        column = config.cmscan_allowed_families_column
        if column in df.columns:
            allowed.update({str(name).strip().lower() for name in df[column].dropna().unique()})
    return allowed


def parse_cmscan_tblout(tblout_path: Path) -> pd.DataFrame:
    if not tblout_path.exists():
        return pd.DataFrame()
    rows: List[Dict[str, Any]] = []
    with open(tblout_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 17:
                continue
            try:
                row = {
                    "target_name": parts[0],
                    "target_accession": parts[1],
                    "sequence_id": parts[2],
                    "sequence_accession": parts[3],
                    "model_type": parts[4],
                    "model_start": int(parts[5]),
                    "model_end": int(parts[6]),
                    "sequence_start": int(parts[7]),
                    "sequence_end": int(parts[8]),
                    "strand": parts[9],
                    "truncation": parts[10],
                    "pass": parts[11],
                    "gc": float(parts[12]),
                    "bias": float(parts[13]),
                    "score": float(parts[14]),
                    "e_value": float(parts[15]),
                    "inc": parts[16],
                    "description": " ".join(parts[17:]) if len(parts) > 17 else "",
                }
            except ValueError:
                continue
            rows.append(row)
    return pd.DataFrame(rows)


def run_cmscan_filter(candidate_table: pd.DataFrame, config: Config) -> pd.DataFrame:
    if candidate_table.empty or not config.cmscan_model_paths:
        empty = pd.DataFrame()
        empty.to_csv(config.cmscan_hits_csv, index=False)
        return empty
    if not config.rna_candidates_fasta.exists():
        empty = pd.DataFrame()
        empty.to_csv(config.cmscan_hits_csv, index=False)
        return empty

    hits_frames: List[pd.DataFrame] = []
    for idx, cm_path in enumerate(config.cmscan_model_paths):
        if not cm_path.exists():
            continue
        tblout_path = config.cmscan_tblout_dir / f"cmscan_{cm_path.stem}_{idx}.tblout"
        cmd = [
            config.cmscan_binary,
            "--tblout",
            str(tblout_path),
            str(cm_path),
            str(config.rna_candidates_fasta),
        ]
        subprocess.run(cmd, check=True)
        tbl_df = parse_cmscan_tblout(tblout_path)
        if not tbl_df.empty:
            tbl_df["cm_source"] = str(cm_path)
            hits_frames.append(tbl_df)

    if not hits_frames:
        return pd.DataFrame()

    hits = pd.concat(hits_frames, ignore_index=True)
    allowed_families = load_allowed_families(config)
    if allowed_families:
        hits = hits[hits["target_name"].str.lower().isin(allowed_families)]

    if hits.empty:
        hits.to_csv(config.cmscan_hits_csv, index=False)
        return hits

    hits = hits[hits["e_value"] <= config.cmscan_evalue_threshold]
    if hits.empty:
        hits.to_csv(config.cmscan_hits_csv, index=False)
        return hits

    seq_to_root = dict(
        zip(candidate_table["Sequence_ID"].astype(str), candidate_table["Root ID"].astype(str))
    )
    hits["Root ID"] = hits["sequence_id"].map(seq_to_root)
    hits = hits.dropna(subset=["Root ID"])
    hits.to_csv(config.cmscan_hits_csv, index=False)
    return hits


def run_pipeline(config_path: Path) -> None:
    config = load_config(config_path)
    prompt_seqs = read_prompts(config.input_prompts, config.batched, config.batch_size)
    model, tokenizer = model_load(config.model_name)
    prompts, sequences, scores, ids = sample_model(
        prompt_batches=prompt_seqs,
        model=model,
        tokenizer=tokenizer,
        file_save_location=str(config.evo_gen_seqs_file_save_location),
        n_tokens=config.n_tokens,
        temp=config.temperature,
        top_k=config.top_k,
        batched=config.batched,
        n_sample_per_prompt=config.n_sample_per_prompt,
        force_prompt_threshold=2,
    )

    final_sequences = get_rc(sequences, rc_truth=config.rc_truth, return_both=config.return_both)
    make_fasta(final_sequences, prompts, ids, str(config.all_seqs_fasta))
    run_prodigal(str(config.all_seqs_fasta), str(config.proteins_file), str(config.orfs_file))

    filter_protein_fasta(
        str(config.proteins_file),
        str(config.filtered_proteins_file),
        str(config.segmasker_path),
        config.filter_min_length,
        config.filter_max_length,
        config.filter_partial_bool,
        config.segmasker_threshold,
    )

    fold_stats = fold_proteins(str(config.filtered_proteins_file), str(config.output_folds_file))
    if config.run_esm_fold:
        filtered_folds = filter_proteins_by_threshold(
            fold_stats,
            str(config.output_filtered_folds),
            config.plddt_threshold,
            config.ptm_threshold,
        )
    else:
        filtered_folds = fold_stats

    if isinstance(filtered_folds, pd.DataFrame) and not filtered_folds.empty:
        filtered_folds = filtered_folds.copy()
        if "Root ID" not in filtered_folds.columns:
            filtered_folds["Root ID"] = filtered_folds["Evo Sequence ID"].astype(str).str.split("_").str[0]

    trf_df = get_tandem_repeats(filtered_folds, config.evo_gen_seqs_file_save_location, config)
    fold_trf_df = fold_trfs(trf_df, config.rna_fold_csv)
    fold_trf_df, trf_root_ids = filter_folded_trfs(trf_df, fold_trf_df, config)
    if trf_root_ids:
        trf_df = trf_df[trf_df["Root ID"].astype(str).isin(trf_root_ids)].copy()
        trf_df.to_csv(config.output_trf_csv, index=False)
        if isinstance(filtered_folds, pd.DataFrame) and not filtered_folds.empty:
            filtered_folds = filtered_folds[filtered_folds["Root ID"].astype(str).isin(trf_root_ids)].copy()
    else:
        empty_columns = list(trf_df.columns) if isinstance(trf_df, pd.DataFrame) else []
        trf_df = pd.DataFrame(columns=empty_columns)
        trf_df.to_csv(config.output_trf_csv, index=False)
        if isinstance(filtered_folds, pd.DataFrame):
            filtered_folds = pd.DataFrame(columns=filtered_folds.columns)

    ta_pairs_df = get_at_pairs(fold_trf_df, filtered_folds, config.ta_pairs_csv)
    visualize_rna_structures(config.rna_structures_reference_csv, fold_trf_df)

    candidate_table = prepare_rna_candidate_table(trf_df, fold_trf_df, config.rna_candidates_csv)
    if not candidate_table.empty:
        write_rna_candidates_fasta(candidate_table, config.rna_candidates_fasta)
        rna_structure_hits = run_rna_structure_filter(config.rna_candidates_csv, config)
        rna_sequence_hits = run_rna_sequence_filter(config.rna_candidates_csv, config)
        cmscan_hits = run_cmscan_filter(candidate_table, config)
    else:
        if config.rna_candidates_fasta.exists():
            config.rna_candidates_fasta.unlink()
        rna_structure_hits = set()
        rna_sequence_hits = set()
        cmscan_hits = pd.DataFrame()

    hmmscan_hits = run_hmmscan_filter(config)

    hmmscan_root_ids: Set[str] = set()
    if not hmmscan_hits.empty:
        hmmscan_root_ids = set(hmmscan_hits["sequence_id"].astype(str).str.split("_").str[0])

    cmscan_root_ids: Set[str] = set()
    if isinstance(cmscan_hits, pd.DataFrame) and not cmscan_hits.empty:
        cmscan_root_ids = set(cmscan_hits["Root ID"].astype(str))

    keep_root_ids = rna_structure_hits.union(rna_sequence_hits).union(hmmscan_root_ids).union(cmscan_root_ids)

    if isinstance(filtered_folds, pd.DataFrame) and not filtered_folds.empty:
        if keep_root_ids:
            final_candidates = filtered_folds[filtered_folds["Root ID"].astype(str).isin(keep_root_ids)].copy()
        else:
            final_candidates = pd.DataFrame(columns=filtered_folds.columns)
        final_candidates.to_csv(config.final_candidates_csv, index=False)

    if not ta_pairs_df.empty:
        if keep_root_ids:
            ta_pairs_filtered = ta_pairs_df[ta_pairs_df["Root ID"].astype(str).isin(keep_root_ids)].copy()
        else:
            ta_pairs_filtered = ta_pairs_df.iloc[0:0].copy()
        ta_pairs_filtered.to_csv(config.ta_pairs_csv, index=False)

    print("Pipeline completed successfully.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Type III toxin-antitoxin sampling pipeline.")
    parser.add_argument("--config", required=True, help="Path to configuration file (YAML).")
    args = parser.parse_args()
    run_pipeline(Path(args.config))
