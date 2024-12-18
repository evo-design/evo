"""
Gene completion eval pipeline using Evo.

Usage: python gene_completion.py --config <config_file_path>
"""

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple, NamedTuple, Union
from dataclasses import dataclass
import logging
import json
import os
import pandas as pd
import tempfile
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

from semantic_mining import read_prompts, model_load, sample_model, get_rc, run_prodigal, filter_protein_fasta

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

    # Input/Output paths
    input_prompts: str  # CSV containing input prompts
    evo_gen_seqs_file_save_location: (
        str  # CSV continaing sequences generated by Evo and corresponding prompts
    )
    all_seqs_fasta: str  # FASTA containing Evo generated sequences used as input to Prodigal
    proteins_file: str  # FASTA containing all prodigal-detected protein sequences
    orfs_file: str  # FASTA containing all prodigal-detected ORFs
    filtered_proteins_file: str  # FASTA containing all protein sequences post-complexity and length filtering
    reference_seqs: str  # FASTA containing all reference protein sequences if aligment against reference sequences is desired
    msa_filtered_proteins_fasta: (
        str  # FASTA containing protein sequences that pass aligment against reference sequences
    )
    sequence_alignment_csv: str  # CSV containing all alignments for alignment against references
    output_msa_csv: Path  # CSV containing summarized alignments for alignment against references
    output_summary_csv: (
        Path  # CSV containing summary statistics about overall sequence identities for all prompts
    )
    segmasker_path: Path  # Path to Segmasker installation
    mafft_path: Path  # Path to MAFFT installation

    # Model parameters
    model_name: str  # Model name to use for sampling
    n_tokens: int  # Number of tokens of sequence to be generated
    temperature: float  # Temperature to sample at
    top_k: int  # Top k value for sampling
    batched: bool  # If generations should be batched (recommend setting to true for quicker generation)
    batch_size: int  # Size of batch for batched generation
    n_sample_per_prompt: int  # Number of times a sequence should be generated for a given prompt (e.g., 3 means sample 3 times on a given prompt)

    # Filter parameters
    rc_truth: bool  # True if the reverse complement of the generated sequence should be included when downstream processing
    return_both: bool  # True if both the reverse complement and original generated sequence should be used for downstream processing
    filter_min_length: int  # Minimum length of protein to keep during filtering
    filter_max_length: int  # Maximum length of protein to keep during filtering
    filter_partial_bool: bool  # True if partial ORFs should be removed during downstream filtering, only set to True if -p meta flag is not used
    segmasker_threshold: float  # Proportion of protein sequence that can be low-complexity for sequence to be kept during filtering
    run_msa: bool  # True if MSA against a list of reference sequences should be run during filtering
    seq_identity_match_threshold: float  # If running MSA, minimum sequence identity that must exist between generated protein sequence and any reference sequence


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
    """
    dna_seq_records = [
        SeqRecord(Seq(prompt + dna_seq), id=seq_id, description=prompt)
        for dna_seq, seq_id, prompt in zip(sequences, ids, prompts)
    ]
    with open(output_file, "w") as output_handle:
        SeqIO.write(dna_seq_records, output_handle, "fasta")


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

    Generated Files:
        output_csv: CSV file with columns:
            - query_id: ID of query sequence
            - reference_id: ID of best matching reference sequence
            - identity: Percent identity between query and reference

        filtered_fasta: FASTA file containing query sequences that had matches
            above the identity threshold
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
    """
    Calculate sequence identity for the non-prompt portion using MAFFT alignment.
    """
    prompt_dna_trimmed = prompt_dna[: -(len(prompt_dna) % 3)] if len(prompt_dna) % 3 != 0 else prompt_dna
    prompt_aa = str(Seq(prompt_dna_trimmed).translate())

    if prompt_aa not in input_aa or prompt_aa not in reference_aa:
        return 0.0

    input_start = input_aa.index(prompt_aa) + len(prompt_aa)
    ref_start = reference_aa.index(prompt_aa) + len(prompt_aa)

    input_non_prompt = input_aa[input_start:]
    ref_non_prompt = reference_aa[ref_start:]

    identity = calculate_sequence_identity(input_non_prompt, ref_non_prompt, mafft_path)

    return identity


def create_summary_statistics(results_df: pd.DataFrame, output_path: Path) -> None:
    """
    Create summary statistics from gene completion analysis results.
    Args:
        results_df: DataFrame containing sequence analysis results
        output_path: Path to save summary statistics CSV

    Generated Files:
        output_path: CSV file containing grouped summary statistics
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
    """
    Args:
        input_fasta: FASTA file containing generated gene sequences with UUID headers
        uuid_prompts_csv: CSV mapping UUIDs to prompts with columns:
            - UUID: Unique identifier for each generated sequence
            - Prompt: The prompt used to generate sequence (truncated gene)
        prompt_info_csv: CSV with prompt metadata containing:
            - Shortened_Sequence: The truncated gene sequence used as prompt
            - Protein_Label: Label identifying the protein type
            - Length_Percentage: What percentage of gene was provided
        reference_fasta: FASTA file containing full reference gene sequences
        output_csv: Path to save detailed analysis results
        output_summary_csv: Path to save summary statistics

    Generated Files:
        output_csv: CSV file containing sequence analysis with columns:
            - UUID: Unique sequence identifier
            - Input_Sequence: Generated sequence
            - Prompt: Truncated sequence used as prompt
            - Protein_Label: Type/label of the protein
            - Length_Percentage: Percent of gene provided in prompt
            - Reference_Sequence: Full reference sequence
            - Full_Sequence_Identity: Percent identity to reference
            - Non_Prompt_Sequence_Identity: Identity of generated portion
            - Prompt_Length: Length of translated prompt sequence

        output_summary_csv: Summary statistics grouped by Protein_Label and
            Length_Percentage, including means and standard deviations
    """
    input_sequences: Dict[str, str] = {
        record.id.split(" ")[0].split("_")[0]: str(record.seq).replace("*", "")
        for record in SeqIO.parse(input_fasta, "fasta")
    }
    reference_sequences: Dict[str, str] = {
        record.id: str(record.seq) for record in SeqIO.parse(reference_fasta, "fasta")
    }
    uuid_prompts_df = pd.read_csv(uuid_prompts_csv)
    prompt_info_df = pd.read_csv(prompt_info_csv)

    results: List[SequenceResult] = []
    for uuid_val, input_seq in input_sequences.items():
        try:
            prompt_row = uuid_prompts_df[uuid_prompts_df["UUID"] == uuid_val]
            if prompt_row.empty:
                continue

            prompt = prompt_row["Prompt"].iloc[0]

            info_row = prompt_info_df[prompt_info_df["Shortened_Sequence"] == prompt]

            if info_row.empty:
                continue

            result = SequenceResult(
                UUID=uuid_val,
                Input_Sequence=input_seq,
                Prompt=prompt,
                Protein_Label=info_row["Protein_Label"].iloc[0],
                Length_Percentage=info_row["Length_Percentage"].iloc[0],
                Reference_Sequence=reference_sequences.get(info_row["Protein_Label"].iloc[0], ""),
                Full_Sequence_Identity=calculate_sequence_identity(
                    input_seq, reference_sequences.get(info_row["Protein_Label"].iloc[0], ""), mafft_path
                ),
                Non_Prompt_Sequence_Identity=calculate_non_prompt_sequence_identity(
                    input_seq,
                    reference_sequences.get(info_row["Protein_Label"].iloc[0], ""),
                    prompt,
                    mafft_path,
                ),
                Prompt_Length=len(translate_dna_sequence(prompt)),
            )
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
        config_path: Path to JSON configuration file

    Generated Files:
        Multiple output files specified in config:
        - Generated sequences FASTA
        - Protein predictions from Prodigal
        - Filtered protein sequences
        - MSA results
        - Analysis CSV
        - Summary statistics
    """
    # Load configuration
    with open(config_path) as f:
        config_dict = json.load(f)
    config = Config(**config_dict)

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
    filter_protein_fasta(
        config.proteins_file,
        config.filtered_proteins_file,
        config.segmasker_path,
        config.filter_min_length,
        config.filter_max_length,
        config.filter_partial_bool,
        config.segmasker_threshold,
    )
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run sampling script with a configuration file.")
    parser.add_argument(
        "--config", required=True, help="Path to the configuration file (e.g., path/to/config.json)"
    )
    args = parser.parse_args()
    run_pipeline(args.config)