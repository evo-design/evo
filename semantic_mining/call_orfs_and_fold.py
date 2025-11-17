"""
Acr sampling pipeline using Evo.

Usage: python acr_sample.py --config <config_file_path>
"""

import argparse
import json
import pandas as pd
from typing import List, Tuple, Dict, Any, Union, TypeVar
from dataclasses import dataclass
from transformers import (
    AutoTokenizer,
    EsmForProteinFolding,
)
from stripedhyena.model import StripedHyena
from stripedhyena.tokenizer import CharLevelTokenizer

from semantic_mining import (
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

BatchType = List[List[Tuple[str, str, str]]]
PromptType = Union[str, List[str]]
SeqRecord = TypeVar("SeqRecord")
import uuid


@dataclass
class Config:
    """Configuration for protein sequence generation and analysis."""
    # Input/Output paths
    all_seqs_fasta: str  # FASTA containing Evo generated sequences used as input to Prodigal
    proteins_file: str  # FASTA containing all prodigal-detected protein sequences
    orfs_file: str  # FASTA containing all prodigal-detected ORFs
    filtered_proteins_file: str  # FASTA containing all protein sequences post-complexity and length filtering
    segmasker_path: str  # Path to Segmasker installation
    output_folds_file: (
        str  # CSV containing folds and sequence information for all protein sequences folded using ESMFold
    )
    output_filtered_folds: str  # FASTA containing sequences that passed fold quality filtering
    
    # Filter parameters
    rc_truth: bool  # True if the reverse complement of the generated sequence should be included when downstream processing
    return_both: bool  # True if both the reverse complement and original generated sequence should be used for downstream processing
    filter_min_length: int  # Minimum length of protein to keep during filtering
    filter_max_length: int  # Maximum length of protein to keep during filtering
    filter_partial_bool: bool  # True if partial ORFs should be removed during downstream filtering, only set to True if -p meta flag is not used
    segmasker_threshold: float  # Proportion of protein sequence that can be low-complexity for sequence to be kept during filtering
    run_esm_fold: bool  # True if protein sequences should be folded during filtering
    plddt_threshold: float  # Minimum pLDDT value of folded protein sequence to keep during filtering
    ptm_threshold: float  # Minimum pTM value of folded protein sequence to keep during filtering

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> "Config":
        """Create Config instance from dictionary."""
        return cls(**config_dict)


def load_config(config_file: str) -> Config:
    """Load configuration from JSON file."""
    try:
        with open(config_file, "r") as f:
            config_dict = json.load(f)
        return Config.from_dict(config_dict)
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    except json.JSONDecodeError:
        raise ValueError(f"Invalid JSON in configuration file: {config_file}")

def process_sequence_data(csv_path):
    """
    Process sequence data from a CSV file, cleaning sequences and extracting model information.
    
    Parameters:
    csv_path (str): Path to the CSV file
    
    Returns:
    tuple: (prompts, sequences, ids) where each is a list
    """
    # Read the CSV file
    df = pd.read_csv(csv_path, header=None)
    
    # Initialize output lists
    prompts = []
    sequences = []
    ids = []
    
    # Process each row
    for _, row in df.iterrows():
        # Extract sequence (first column)
        sequence = str(row[0])
        # Clean sequence to keep only ATGC characters
        cleaned_sequence = ''.join(char for char in sequence.upper() if char in 'ATGC')
        # Generate UUID
        sequence_id = str(uuid.uuid4())
        
        # Create prompt
        prompt = 'evo_7b_1m_gtdb_m_genitalium'
        
        # Append to output lists
        prompts.append(prompt)
        sequences.append(cleaned_sequence)
        ids.append(sequence_id)
    
    return prompts, sequences, ids

def process_sequences(config: Config) -> None:
    """Process sequences through the pipeline."""
    print("Starting sequence processing...")
    print("Prompts loaded")

    prompts,sequences,ids = process_sequence_data('/large_storage/hielab/userspace/adititm/evo2/misc/m_genitalium.csv')
    # Generate sequences
    # Process sequences
    final_sequences = get_rc(sequences, rc_truth=config.rc_truth, return_both=config.return_both)

    # Create FASTA files
    make_fasta(final_sequences, prompts, ids, config.all_seqs_fasta)

    # Run protein analysis pipeline
    run_prodigal(config.all_seqs_fasta, config.proteins_file, config.orfs_file)

    print("Base protein filtering started...", flush=True)
    filter_protein_fasta(
        config.proteins_file,
        config.filtered_proteins_file,
        config.segmasker_path,
        config.filter_min_length,
        config.filter_max_length,
        config.filter_partial_bool,
        config.segmasker_threshold,
    )
    print("Base protein filtering complete", flush=True)


def process_folds(config: Config) -> pd.DataFrame:
    """Process sequence alignments and protein folding."""
    print("Starting protein folding...", flush=True)
    fold_stats = fold_proteins(config.filtered_proteins_file, config.output_folds_file)
    print("Protein folding complete", flush=True)
    filtered_folds = filter_proteins_by_threshold(
        fold_stats, config.output_filtered_folds, config.plddt_threshold, config.ptm_threshold
    )

    return filtered_folds


def main(config_file: str) -> None:
    try:
        # # Load configuration
        config = load_config(config_file)
        # print("Configuration loaded", flush=True)

        # # Execute pipeline
        # process_sequences(config)
        if config.run_esm_fold:
            filtered_folds = process_folds(config)

        print("Pipeline completed successfully", flush=True)

    except Exception as e:
        print(f"Error in pipeline execution: {str(e)}", flush=True)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run sampling script with a configuration file.")
    parser.add_argument(
        "--config", required=True, help="Path to the configuration file (e.g., path/to/config.json)"
    )
    args = parser.parse_args()
    main(args.config)
