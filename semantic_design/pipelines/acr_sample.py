"""
Acr sampling pipeline using Evo.

Usage: python pipelines/acr_sample.py --config <config_file_path>
"""

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

import pandas as pd
import yaml
from transformers import (
    AutoTokenizer,
    EsmForProteinFolding,
)
from stripedhyena.model import StripedHyena
from evo.tokenizer import CharLevelTokenizer

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
    """Configuration for protein sequence generation and analysis."""

    # Input/Output paths
    input_prompts: Path  # CSV containing input prompts
    output_dir: Path  # Directory where all outputs will be materialized
    segmasker_path: Path  # Path to Segmasker installation

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
    filter_partial_bool: (
        bool  # True if partial ORFs should be removed during downstream filtering
    )
    segmasker_threshold: float  # Proportion of protein sequence that can be low-complexity for sequence to be kept during filtering
    run_esm_fold: bool  # True if protein sequences should be folded during filtering
    plddt_threshold: (
        float  # Minimum pLDDT value of folded protein sequence to keep during filtering
    )
    ptm_threshold: (
        float  # Minimum pTM value of folded protein sequence to keep during filtering
    )

    evo_gen_seqs_file_save_location: Path = field(init=False)
    all_seqs_fasta: Path = field(init=False)
    proteins_file: Path = field(init=False)
    orfs_file: Path = field(init=False)
    filtered_proteins_file: Path = field(init=False)
    output_folds_file: Path = field(init=False)
    output_filtered_folds: Path = field(init=False)

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> "Config":
        """Create Config instance from dictionary."""
        return cls(**config_dict)

    def __post_init__(self) -> None:
        """Resolve user paths and define all derived output file locations."""
        self.input_prompts = Path(self.input_prompts)
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.segmasker_path = Path(self.segmasker_path)

        self.evo_gen_seqs_file_save_location = (
            self.output_dir / "generated_sequences.csv"
        )
        self.all_seqs_fasta = self.output_dir / "all_sequences.fasta"
        self.proteins_file = self.output_dir / "proteins.fasta"
        self.orfs_file = self.output_dir / "orfs.fasta"
        self.filtered_proteins_file = self.output_dir / "filtered_proteins.fasta"
        self.output_folds_file = self.output_dir / "folds.csv"
        self.output_filtered_folds = self.output_dir / "filtered_folds.csv"


def load_config(config_file: str) -> Config:
    """Load configuration from a YAML file path.

    Args:
        config_file: Path to the user-provided YAML config.

    Returns:
        Config instance with derived output paths created.
    """
    try:
        with open(config_file, "r") as f:
            config_dict = yaml.safe_load(f)
        if not isinstance(config_dict, dict):
            raise ValueError(
                f"Configuration file must contain a YAML mapping: {config_file}"
            )
        return Config.from_dict(config_dict)
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML in configuration file: {config_file}") from exc


def process_sequences(
    config: Config, model: StripedHyena, tokenizer: CharLevelTokenizer
) -> None:
    """Sample DNA sequences with Evo and run the baseline filtering cascade.

    Args:
        config: Loaded configuration for this run.
        model: Language model used for sequence generation.
        tokenizer: Tokenizer compatible with the provided model.

    Returns:
        None. Artifacts are written to the configured output directory.
    """
    print("Starting sequence processing...")

    # Read prompts
    prompt_seqs = read_prompts(config.input_prompts, config.batched, config.batch_size)
    print("Prompts loaded")

    # Generate sequences
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
    # Process sequences
    final_sequences = get_rc(
        sequences, rc_truth=config.rc_truth, return_both=config.return_both
    )

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
    """Fold filtered proteins with ESMFold and return sequences passing thresholds.

    Args:
        config: Run configuration containing path and cutoff information.

    Returns:
        DataFrame of proteins that met the user-defined pLDDT/pTM thresholds.
    """
    print("Starting protein folding...", flush=True)
    fold_stats = fold_proteins(config.filtered_proteins_file, config.output_folds_file)
    print("Protein folding complete", flush=True)
    filtered_folds = filter_proteins_by_threshold(
        fold_stats,
        config.output_filtered_folds,
        config.plddt_threshold,
        config.ptm_threshold,
    )

    return filtered_folds


def main(config_file: str) -> None:
    """Entry point for running the Acr sampling pipeline via a config file path.

    Args:
        config_file: Path to the YAML configuration.

    Returns:
        None.
    """
    try:
        # Load configuration
        config = load_config(config_file)
        print("Configuration loaded", flush=True)

        model, tokenizer = model_load(config.model_name)

        # Execute pipeline
        process_sequences(config, model, tokenizer)
        if config.run_esm_fold:
            filtered_folds = process_folds(config)

        print("Pipeline completed successfully", flush=True)

    except Exception as e:
        print(f"Error in pipeline execution: {str(e)}", flush=True)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run sampling script with a configuration file."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to the configuration file (e.g., path/to/config.yaml)",
    )
    args = parser.parse_args()
    main(args.config)
