"""
Gene completion eval pipeline using Evo.

Usage: python pipelines/gene_completion.py --config <config_file_path>
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple, Union, Set


import pandas as pd
import yaml
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from semantic_design import read_prompts, model_load, sample_model, get_rc, run_prodigal, filter_protein_fasta

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


# Type definitions
class ModelOutput(NamedTuple):
    sequences: List[str]
    scores: List[float]


class BatchData(NamedTuple):
    prompts: List[str]
    sequences: List[str]
    scores: List[str]
    ids: List[str]


@dataclass
class SequenceResult:
    """Data class to store sequence analysis results."""

    UUID: str
    Input_Sequence: str
    Prompt: str
    Protein_Label: str
    Length_Percentage: float
    Reference_Sequence: str
    Full_Sequence_Identity: Optional[float]
    Non_Prompt_Sequence_Identity: Optional[float]
    Prompt_Length: int


@dataclass
class Config:
    """Configuration for protein sequence generation and analysis."""

    # Input/Output paths supplied by user
    input_prompts: Path  # CSV containing input prompts
    reference_seqs: Path  # FASTA of reference protein sequences for optional alignment
    output_dir: Path  # Directory that will contain all generated artifacts
    segmasker_path: Path  # Path to Segmasker installation
    mafft_path: Path  # Path to MAFFT installation

    # Model parameters
    model_name: str
    n_tokens: int
    temperature: float
    top_k: int
    batched: bool
    batch_size: int
    n_sample_per_prompt: int

    # Filter parameters
    rc_truth: bool
    return_both: bool
    filter_min_length: int
    filter_max_length: int
    filter_partial_bool: bool
    segmasker_threshold: float
    run_msa: bool
    seq_identity_match_threshold: float

    # Derived output paths
    evo_gen_seqs_file_save_location: Path = field(init=False)
    all_seqs_fasta: Path = field(init=False)
    proteins_file: Path = field(init=False)
    orfs_file: Path = field(init=False)
    filtered_proteins_file: Path = field(init=False)
    msa_filtered_proteins_fasta: Path = field(init=False)
    sequence_alignment_csv: Path = field(init=False)
    output_msa_csv: Path = field(init=False)
    output_summary_csv: Path = field(init=False)

    def __post_init__(self) -> None:
        """Normalize input paths and establish every derived output artifact location."""
        self.input_prompts = Path(self.input_prompts)
        self.reference_seqs = Path(self.reference_seqs)
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.segmasker_path = Path(self.segmasker_path)
        self.mafft_path = Path(self.mafft_path)

        self.evo_gen_seqs_file_save_location = self.output_dir / "generated_sequences.csv"
        self.all_seqs_fasta = self.output_dir / "all_sequences.fasta"
        self.proteins_file = self.output_dir / "proteins.fasta"
        self.orfs_file = self.output_dir / "orfs.fasta"
        self.filtered_proteins_file = self.output_dir / "filtered_proteins.fasta"
        self.msa_filtered_proteins_fasta = self.output_dir / "msa_filtered_proteins.fasta"
        self.sequence_alignment_csv = self.output_dir / "sequence_alignment.csv"
        self.output_msa_csv = self.output_dir / "msa_results.csv"
        self.output_summary_csv = self.output_dir / "summary_statistics.csv"


def load_config(config_path: Path) -> Config:
    """Load a YAML configuration file into a Config object.

    Args:
        config_path: Path to the YAML config.

    Returns:
        Parsed Config instance with derived paths created.
    """
    config_path = Path(config_path)
    try:
        with open(config_path, "r") as handle:
            data = yaml.safe_load(handle)
        if not isinstance(data, dict):
            raise ValueError(f"Configuration must be a mapping: {config_path}")
        return Config(**data)
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML in configuration file: {config_path}") from exc


def translate_dna_sequence(dna_seq: str) -> str:
    """
     Translate a DNA sequence to protein sequence.
        Args:
        dna_seq: DNA sequence string

    Returns:
        Protein sequence
    """
    trimmed_length = len(dna_seq) - (len(dna_seq) % 3)
    trimmed_seq = dna_seq[:trimmed_length]
    return str(Seq(trimmed_seq).translate())


def make_gene_completion_fasta(
    sequences: List[str], prompts: List[str], ids: List[str], output_file: Path
) -> None:
    """
     Creates a FASTA file with prompts appemded to generations.

    Args:
        sequences: List of DNA sequences
        prompts: List of prompts used to generate the sequences
        ids: List of unique identifiers for each sequence
        output_file: Path to save the output FASTA file

    Generated Files:
        output_file: FASTA format file where each record contains:
            - Header: >sequence_id prompt_text
            - Sequence: prompt_sequence + dna_sequence

    Returns:
        None. Contents are written to `output_file`.
    """
    dna_seq_records = [
        SeqRecord(Seq(prompt + dna_seq), id=seq_id, description=prompt)
        for dna_seq, seq_id, prompt in zip(sequences, ids, prompts)
    ]
    with open(output_file, "w") as output_handle:
        SeqIO.write(dna_seq_records, output_handle, "fasta")


def filter_orfs_by_prompt(
    proteins_fasta: Path, orfs_fasta: Path, prompts_csv: Path
) -> None:
    """Keep only ORFs (and their proteins) whose nucleotide sequence contains the prompt DNA.

    Args:
        proteins_fasta: Path to Prodigal-translated protein FASTA.
        orfs_fasta: Path to Prodigal nucleotide ORF FASTA.
        prompts_csv: CSV with UUID/prompt pairs (`generated_sequences.csv`).

    Returns:
        None. The provided FASTA files are overwritten with filtered content.
    """
    if not Path(orfs_fasta).exists() or not Path(proteins_fasta).exists():
        logger.warning("ORF/protein FASTA missing; skipping prompt-based ORF filtering.")
        return
    if not Path(prompts_csv).exists():
        logger.warning("Prompts CSV missing; skipping prompt-based ORF filtering.")
        return

    prompts_df = pd.read_csv(prompts_csv)
    if "UUID" not in prompts_df.columns or "Prompt" not in prompts_df.columns:
        logger.warning("Prompts CSV missing UUID/Prompt columns; skipping filter.")
        return

    prompt_map = {
        str(row["UUID"]).split("_")[0]: str(row["Prompt"]).upper()
        for _, row in prompts_df.iterrows()
    }

    allowed_ids: Set[str] = set()
    filtered_orfs: List[SeqRecord] = []

    for record in SeqIO.parse(str(orfs_fasta), "fasta"):
        base_id = record.id.split(" ")[0]
        uuid = base_id.split("_")[0]
        prompt_dna = prompt_map.get(uuid)
        if not prompt_dna:
            continue
        seq = str(record.seq).upper()
        if prompt_dna not in seq:
            continue
        if base_id in allowed_ids:
            continue
        allowed_ids.add(base_id)
        filtered_orfs.append(record)

    if not filtered_orfs:
        logger.warning("No ORFs contained their prompts; downstream outputs will be empty.")

    SeqIO.write(filtered_orfs, str(orfs_fasta), "fasta")

    filtered_proteins = [
        record
        for record in SeqIO.parse(str(proteins_fasta), "fasta")
        if record.id.split(" ")[0] in allowed_ids
    ]
    SeqIO.write(filtered_proteins, str(proteins_fasta), "fasta")
    logger.info("Retained %d ORFs/proteins whose nucleotide sequence contains the prompt.", len(filtered_orfs))


def build_reference_lookup(reference_fasta: Path) -> Dict[str, str]:
    """Create a case-insensitive mapping from reference labels to sequences.

    Args:
        reference_fasta: FASTA file containing reference protein sequences.

    Returns:
        Dictionary keyed by various identifiers (accession, description tokens) to amino-acid sequences.
    """
    lookup: Dict[str, str] = {}
    for record in SeqIO.parse(reference_fasta, "fasta"):
        seq = str(record.seq)
        description = record.description.lower()
        candidates = {
            record.id.lower(),
            description,
        }
        candidates.update(token.strip("[](),") for token in description.replace("/", " ").split())
        for key in candidates:
            if key and key not in lookup:
                lookup[key] = seq
    return lookup


def calculate_sequence_identity(seq1: str, seq2: str, mafft_path: str = "mafft") -> Optional[float]:
    """
    Calculate sequence identity between two sequences using MAFFT.

    Args:
        seq1: First amino acid sequence
        seq2: Second amino acid sequence
        mafft_path: Path to MAFFT executable (default: "mafft")

    Returns:
        Percent identity between sequences (0-100), or None if alignment fails
    """
    if not seq1 or not seq2:
        return None

    record1 = SeqRecord(Seq(seq1), id="seq1")
    record2 = SeqRecord(Seq(seq2), id="seq2")

    try:
        aligned_seq1, aligned_seq2, identity = align_pair(record1, record2, mafft_path)
        return identity * 100

    except subprocess.CalledProcessError as e:
        print(f"MAFFT alignment failed: {str(e)}")
        return None


def align_pair(query_record: SeqRecord, ref_record: SeqRecord, mafft_path: str) -> Tuple[str, str, float]:
    """
    Align a pair of sequences using MAFFT.

    Args:
        query_record: First sequence as SeqRecord
        ref_record: Second sequence as SeqRecord
        mafft_path: Path to MAFFT executable

    Returns:
        Tuple containing:
            - Aligned first sequence (string)
            - Aligned second sequence (string)
            - Sequence identity (float between 0-1)
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp_fasta:
        SeqIO.write([query_record, ref_record], tmp_fasta, "fasta")
        tmp_fasta_name = tmp_fasta.name

    try:
        result = subprocess.run([mafft_path, tmp_fasta_name], capture_output=True, text=True, check=True)
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as aligned_file:
            aligned_file.write(result.stdout)
            aligned_file_name = aligned_file.name

        alignment = AlignIO.read(aligned_file_name, "fasta")
        aligned_seq1, aligned_seq2 = alignment[0].seq, alignment[1].seq

        identity_count = sum(
            1 for a, b in zip(aligned_seq1, aligned_seq2) if a != "-" and b != "-" and a == b
        )
        aligned_length = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != "-" and b != "-")
        identity = identity_count / aligned_length if aligned_length > 0 else 0

        return str(aligned_seq1), str(aligned_seq2), identity

    finally:
        for fname in [tmp_fasta_name, aligned_file_name]:
            try:
                os.remove(fname)
            except (OSError, UnboundLocalError):
                pass


def align_and_save_closest_match(
    input_fasta: Path,
    reference_fasta: Path,
    output_csv: Path,
    filtered_fasta: Path,
    identity_threshold: float,
    mafft_path: str = "mafft",
) -> None:
    """
    Align sequences and save closest matches using MAFFT.

    Args:
        input_fasta: Path to input FASTA file containing query sequences
        reference_fasta: Path to reference FASTA file containing sequences to match against
        output_csv: Path to save alignment results CSV
        filtered_fasta: Path to save filtered sequences FASTA
        identity_threshold: Minimum percent identity threshold for keeping matches (0-100)
        mafft_path: Path to MAFFT executable (default: "mafft")

)    Generated Files:
        output_csv: CSV file with columns:
            - query_id: ID of query sequence
            - reference_id: ID of best matching reference sequence
            - identity: Percent identity between query and reference

        filtered_fasta: FASTA file containing query sequences that had matches
            above the identity threshold

    Returns:
        None. Artifacts are written to `output_csv` and `filtered_fasta`.
    """
    reference_seqs = {record.id: record for record in SeqIO.parse(reference_fasta, "fasta")}

    results = []
    filtered_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        best_identity = 0.0
        best_match = None

        for ref_id, ref_record in reference_seqs.items():
            _, _, identity = align_pair(record, ref_record, mafft_path)
            identity *= 100
            if identity > best_identity:
                best_identity = identity
                best_match = ref_id

        if best_identity >= identity_threshold:
            results.append({"query_id": record.id, "reference_id": best_match, "identity": best_identity})
            filtered_records.append(record)

    pd.DataFrame(results).to_csv(output_csv, index=False)
    SeqIO.write(filtered_records, filtered_fasta, "fasta")


def calculate_sequence_identity_w_prompt(seq1: str, seq2: str, mafft_path: str = "mafft") -> float:
    """
    Calculate sequence identity between two sequences using MAFFT.

    Args:
        seq1: First amino acid sequence
        seq2: Second amino acid sequence
        mafft_path: Path to MAFFT executable

    Returns:
        Percent identity between sequences (0-100), defaulting to 0 if alignment fails
    """
    if not seq1 or not seq2:
        return 0.0

    tmp_fasta = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as f:
            tmp_fasta = f.name
            SeqIO.write([SeqRecord(Seq(seq1), id="seq1"), SeqRecord(Seq(seq2), id="seq2")], f, "fasta")

        try:
            result = subprocess.run(
                [mafft_path, "--quiet", tmp_fasta], capture_output=True, text=True, check=True
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            return 0.0

        aligned = list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))
        if len(aligned) != 2:
            return 0.0

        aligned_seq1 = str(aligned[0].seq)
        aligned_seq2 = str(aligned[1].seq)

        matches = sum(c1 == c2 for c1, c2 in zip(aligned_seq1, aligned_seq2))
        return (matches / len(aligned_seq1)) * 100

    finally:
        if tmp_fasta and os.path.exists(tmp_fasta):
            try:
                os.unlink(tmp_fasta)
            except OSError:
                pass


def calculate_non_prompt_sequence_identity(
    input_aa: str, reference_aa: str, prompt_dna: str, mafft_path: str = "mafft"
) -> float:
    """Calculate sequence identity for the non-prompt portion using aligned sequences.

    Args:
        input_aa: Generated full-length amino-acid sequence.
        reference_aa: Reference amino-acid sequence.
        prompt_dna: DNA string representing the prompt region.
        mafft_path: Path to MAFFT (default: "mafft").

    Returns:
        Percent identity (0-100) for the aligned region after the prompt.
    """
    if not input_aa or not reference_aa or not prompt_dna:
        return 0.0

    prompt_dna_trimmed = prompt_dna[: -(len(prompt_dna) % 3)] if len(prompt_dna) % 3 != 0 else prompt_dna
    prompt_aa = str(Seq(prompt_dna_trimmed).translate())
    prompt_len = len(prompt_aa)
    if prompt_len == 0:
        return 0.0

    try:
        aligned_seq1, aligned_seq2, _ = align_pair(
            SeqRecord(Seq(input_aa), id="input"),
            SeqRecord(Seq(reference_aa), id="reference"),
            mafft_path,
        )
    except subprocess.CalledProcessError:
        return 0.0

    consumed_input = 0
    consumed_ref = 0
    matches = 0
    positions = 0

    for aa1, aa2 in zip(aligned_seq1, aligned_seq2):
        if aa1 != "-":
            consumed_input += 1
        if aa2 != "-":
            consumed_ref += 1

        if consumed_input <= prompt_len or consumed_ref <= prompt_len:
            continue

        if aa1 == "-" or aa2 == "-":
            continue

        positions += 1
        if aa1 == aa2:
            matches += 1

    if positions == 0:
        return 0.0
    return (matches / positions) * 100.0


def create_summary_statistics(results_df: pd.DataFrame, output_path: Path) -> None:
    """
    Create summary statistics from gene completion analysis results.
    Args:
        results_df: DataFrame containing sequence analysis results
        output_path: Path to save summary statistics CSV

    Generated Files:
        output_path: CSV file containing grouped summary statistics

    Returns:
        None.
    """
    if results_df.empty:
        logger.error("No results to process. Exiting.")
        return

    results_df["Non_Prompt_Sequence_Identity"] = pd.to_numeric(
        results_df["Non_Prompt_Sequence_Identity"], errors="coerce"
    )

    summary_stats = results_df.groupby(["Prompt", "Protein_Label", "Length_Percentage"]).agg(
        avg_full_identity=("Full_Sequence_Identity", lambda x: x[x > 0].mean()),
        std_full_identity=("Full_Sequence_Identity", lambda x: x[x > 0].std()),
        count_full_identity=("Full_Sequence_Identity", lambda x: x[x > 0].count()),
        avg_non_prompt_identity=("Non_Prompt_Sequence_Identity", lambda x: x[x > 0].mean()),
        std_non_prompt_identity=("Non_Prompt_Sequence_Identity", lambda x: x[x > 0].std()),
        count_non_prompt_identity=("Non_Prompt_Sequence_Identity", lambda x: x[x > 0].count()),
        prompt_length=("Prompt_Length", "first"),
    )

    summary_stats = summary_stats.reset_index()

    summary_stats = summary_stats.fillna(0)

    numeric_cols = [
        "avg_full_identity",
        "std_full_identity",
        "avg_non_prompt_identity",
        "std_non_prompt_identity",
    ]
    summary_stats[numeric_cols] = summary_stats[numeric_cols].round(2)

    summary_stats.to_csv(output_path, index=False)


def process_gene_completion_sequences(
    input_fasta: Path,
    uuid_prompts_csv: Path,
    prompt_info_csv: Path,
    reference_fasta: Path,
    output_csv: Path,
    output_summary_csv: Path,
    mafft_path: Path,
) -> None:
    """Analyze generated sequences against references and write per-gene metrics.

    Args:
        input_fasta: FASTA file containing generated gene sequences with UUID headers.
        uuid_prompts_csv: CSV mapping sequence UUIDs to the prompts used for generation.
        prompt_info_csv: CSV with prompt metadata (prompt string, protein label, length%).
        reference_fasta: FASTA file containing full-length reference gene sequences.
        output_csv: Path to write the per-sequence MSA/identity statistics.
        output_summary_csv: Path to write grouped summary statistics.
        mafft_path: Location of the MAFFT executable.

    Returns:
        None. Results are written to `output_csv` and `output_summary_csv`.
    """
    input_sequences: Dict[str, str] = {
        record.id.split(" ")[0].split("_")[0]: str(record.seq).replace("*", "")
        for record in SeqIO.parse(input_fasta, "fasta")
    }
    reference_lookup = build_reference_lookup(reference_fasta)
    uuid_prompts_df = pd.read_csv(uuid_prompts_csv)
    prompt_info_df = pd.read_csv(prompt_info_csv)

    results: List[SequenceResult] = []
    for uuid_val, input_seq in input_sequences.items():
        try:
            prompt_row = uuid_prompts_df[uuid_prompts_df["UUID"] == uuid_val]
            if prompt_row.empty:
                continue

            prompt = prompt_row["Prompt"].iloc[0]
            prompt_aa = translate_dna_sequence(prompt)

            info_row = prompt_info_df[prompt_info_df["Shortened_Sequence"] == prompt]

            if info_row.empty:
                continue

            if not input_seq.startswith(prompt_aa):
                logger.debug(
                    "Skipping UUID %s because translated prompt does not match sequence start.", uuid_val
                )
                continue

            result = SequenceResult(
                UUID=uuid_val,
                Input_Sequence=input_seq,
                Prompt=prompt,
                Protein_Label=info_row["Protein_Label"].iloc[0],
                Length_Percentage=info_row["Length_Percentage"].iloc[0],
                Reference_Sequence="",
                Full_Sequence_Identity=0.0,
                Non_Prompt_Sequence_Identity=0.0,
                Prompt_Length=len(prompt_aa),
            )

            ref_label = str(info_row["Protein_Label"].iloc[0]).lower()
            reference_seq = reference_lookup.get(ref_label)
            if not reference_seq:
                logger.warning("No reference sequence found for label '%s'; skipping.", ref_label)
                continue

            result.Reference_Sequence = reference_seq

            full_identity = calculate_sequence_identity(input_seq, reference_seq, mafft_path) or 0.0
            non_prompt_identity = calculate_non_prompt_sequence_identity(
                input_seq,
                reference_seq,
                prompt,
                mafft_path,
            )
            non_prompt_identity = non_prompt_identity if non_prompt_identity is not None else 0.0

            result.Full_Sequence_Identity = full_identity
            result.Non_Prompt_Sequence_Identity = non_prompt_identity
            results.append(result)

        except Exception as e:
            logger.error(f"Error processing UUID {uuid_val}: {str(e)}")
            continue

    output_df = pd.DataFrame([vars(r) for r in results])
    output_df.to_csv(output_csv, index=False)
    create_summary_statistics(output_df, output_summary_csv)


def run_pipeline(config_path: Path) -> None:
    """
    Execute full gene completion analysis pipeline.
    Args:
        config_path: Path to YAML configuration file

    Generated Files:
        Multiple output files specified in config:
        - Generated sequences FASTA
        - Protein predictions from Prodigal
        - Filtered protein sequences
        - MSA results
        - Analysis CSV
        - Summary statistics

    Returns:
        None.
    """
    config = load_config(config_path)

    # Generate sequences
    prompt_seqs = read_prompts(config.input_prompts, config.batched, config.batch_size)
    model, tokenizer = model_load(config.model_name)

    batch_data = sample_model(
        prompt_batches=prompt_seqs,
        model=model,
        tokenizer=tokenizer,
        file_save_location=config.evo_gen_seqs_file_save_location,
        n_tokens=config.n_tokens,
        temp=config.temperature,
        top_k=config.top_k,
        batched=config.batched,
        n_sample_per_prompt=config.n_sample_per_prompt,
        force_prompt_threshold=2,
    )
    prompts, sequences, scores, ids = batch_data
    # Generate and save FASTA files
    final_sequences = get_rc(sequences, rc_truth=config.rc_truth, return_both=config.return_both)
    make_gene_completion_fasta(final_sequences, prompts, ids, config.all_seqs_fasta)

    # Process sequences
    run_prodigal(config.all_seqs_fasta, config.proteins_file, config.orfs_file)
    filter_orfs_by_prompt(
        config.proteins_file, config.orfs_file, config.evo_gen_seqs_file_save_location
    )
    filter_protein_fasta(
        config.proteins_file,
        config.filtered_proteins_file,
        config.segmasker_path,
        config.filter_min_length,
        config.filter_max_length,
        config.filter_partial_bool,
        config.segmasker_threshold,
    )
    if config.run_msa:
        align_and_save_closest_match(
            config.filtered_proteins_file,
            config.reference_seqs,
            config.sequence_alignment_csv,
            config.msa_filtered_proteins_fasta,
            config.seq_identity_match_threshold,
            config.mafft_path,
        )
        process_gene_completion_sequences(
            input_fasta=config.msa_filtered_proteins_fasta,
            uuid_prompts_csv=config.evo_gen_seqs_file_save_location,
            prompt_info_csv=config.input_prompts,
            reference_fasta=config.reference_seqs,
            output_csv=config.output_msa_csv,
            output_summary_csv=config.output_summary_csv,
            mafft_path=config.mafft_path,
        )
    else:
        logger.info("Skipping MSA and summary generation because run_msa is set to False.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run sampling script with a configuration file.")
    parser.add_argument(
        "--config", required=True, help="Path to the configuration file (e.g., path/to/config.yaml)"
    )
    args = parser.parse_args()
    run_pipeline(args.config)
