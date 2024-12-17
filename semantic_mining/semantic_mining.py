import csv
import json
import sys
import uuid
import math
import subprocess
from pathlib import Path
from typing import Dict, Any, Union, TypeVar, List, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import subprocess
import tempfile
import os
import torch

from stripedhyena.model import StripedHyena
from stripedhyena.tokenizer import CharLevelTokenizer
from evo import Evo
from evo.generation import generate

from transformers import AutoConfig, AutoModelForCausalLM, AutoTokenizer, EsmForProteinFolding, set_seed

BatchType = List[List[Tuple[str, str, str]]]
PromptType = Union[str, List[str]]


def read_prompts(
    input_file: str, batched: bool = True, batch_size: int = 150
) -> Union[List[List[str]], List[str]]:
    """
    Read protein sequences from a CSV file and optionally batch them by sequence length.

    Args:
        input_file: Path to the input CSV file containing sequences.
                   Expected format: First column contains sequences, first row is header.
                   Example:
                   Sequence
                   ATGCATTTTT...
                   GCCGGAATTA...

        batched: If True, sequences are grouped by length and batched.
                If False, returns a simple list of sequences.
                Default is True for optimal model processing.

        batch_size: Maximum number of sequences in each batch when batched=True.
                   Only applies when batched=True.
                   Default is 150 sequences per batch.

    Returns:
        If batched=True:
            List[List[str]]: List of batches, where each batch contains sequences
                            of the same length. Each inner list has max batch_size sequences.
        If batched=False:
            List[str]: Simple list of all sequences.

    Raises:
        FileNotFoundError: If input_file doesn't exist
        csv.Error: If CSV file is malformed
    """
    # Read all sequences from CSV file
    with open(input_file, encoding="utf-8-sig", newline="") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        promptseqs = [row[0] for row in reader]

    # If batching not required, return simple list of sequences
    if not batched:
        return promptseqs

    # Dictionary to group sequences by their length
    length_to_seqs: Dict[int, List[str]] = {}
    for seq in promptseqs:
        seq_length = len(seq)
        if seq_length not in length_to_seqs:
            length_to_seqs[seq_length] = []
        length_to_seqs[seq_length].append(seq)

    # Create batches of sequences with same length
    batches = []
    for seq_data in length_to_seqs.values():
        if len(seq_data) > 1:
            for i in range(0, len(seq_data), batch_size):
                batch = seq_data[i : i + batch_size]
                batches.append(batch)
        else:
            batches.append(seq_data)

    return batches


def model_load(model_name: str) -> tuple[StripedHyena, CharLevelTokenizer]:
    """
    Load model based on provided model name

    Args:
        input_file: Name of model to use for sampling

    Returns:
        A tuple containing:
            - Model loaded using StripedHyena
            - Tokenizer
    """
    evo_model = Evo(model_name, device="cuda")
    evo_model.model = evo_model.model.to("cuda")
    model, tokenizer = evo_model.model, evo_model.tokenizer
    return model, tokenizer


def run_model(
    prompts: PromptType,
    model: StripedHyena,
    tokenizer: CharLevelTokenizer,
    n_tokens: int = 1000,
    temp: float = 0.7,
    top_k: int = 4,
    batched: bool = True,
    device: str = "cuda:0",
    force_prompt_threshold: int = 2,
) -> Tuple[Union[List[str], str], Union[List[float], float]]:
    """
    Generate DNA sequences using Evo.

    Args:
        prompts: Input DNA sequence(s) to seed the generation.
                Can be either a single string or a list of strings when batched=True

        model: Initialized Evo model for sequence generation

        tokenizer: CharLevelTokenizer for converting sequences to/from tokens

        n_tokens: Number of tokens to generate for each sequence

        temp: Temperature parameter for controlling randomness in generation
              Higher values (e.g. 1.0) = more diverse but potentially less reliable sequences
              Lower values (e.g. 0.1) = more conservative but potentially repetitive sequences

        top_k: Number of highest probability tokens to sample from at each step
              Higher values allow more diversity, lower values are more conservative

        batched: Whether to process prompts as a batch (True) or individually (False)
                Batching is more efficient but requires prompts of same length

        device: Device to run generation on ('cuda:0', 'cuda:1', 'cpu', etc)

        force_prompt_threshold: Minimum number of tokens that must match the prompt

    Returns:
        A tuple containing:
        - Generated sequences (str if single sequence, List[str] if batched)
        - Generation scores (float if single sequence, List[float] if batched)
    """
    return generate(
        prompt_seqs=prompts,
        model=model,
        tokenizer=tokenizer,
        n_tokens=n_tokens,
        temperature=temp,
        top_k=top_k,
        batched=batched,
        device=device,
        force_prompt_threshold=force_prompt_threshold,
        cached_generation=True,
        verbose=True,
    )


def read_evo_seqs(gen_seqs: List[str]) -> Tuple[List[str], List[str], List[float], List[str]]:
    """
    Read in Evo-generated sequences.

    Args:
       gen_seqs: List of sequences generated by Evo where each element contains:
                [prompt, generated_sequence, score, sequence_id]
    Returns:
       A tuple containing four parallel lists:
       - prompts (List[str]): Original input sequences used as prompts
       - sequences (List[str]): Generated DNA sequences
       - scores (List[float]): Generation scores (converted to float)
       - ids (List[str]): Unique identifiers for each sequence

    """
    scores = []
    for row in gen_seqs:
        try:
            scores.append(float(row[3]))
        except (ValueError, TypeError):
            scores.append(0.0)

    return (
        [row[1] for row in gen_seqs],  # prompts
        [row[2] for row in gen_seqs],  # sequences
        scores,  # converted scores
        [row[0] for row in gen_seqs],  # UUIDs
    )


def get_rc(sequences: List[str], rc_truth: bool = True, return_both: bool = True) -> List[Seq]:
    """
    Generate reverse complements of DNA sequences

    Args:
        sequences: List of DNA sequences as strings (e.g., ["ATCG", "GCTA"])
        rc_truth: If True, generate reverse complements. If False, only return
                 original sequences.
        return_both: If True and rc_truth is True, return both original sequences
                    and their reverse complements. If False and rc_truth is True,
                    return only reverse complements.
    Returns:
        List[Seq]: A list of Biopython Seq objects containing:
        - If rc_truth=True and return_both=True:
          [reverse_complements + original_sequences]
        - If rc_truth=True and return_both=False:
          [reverse_complements]
        - If rc_truth=False:
          [original_sequences]

    """
    dna_seq = [Seq(curr_seq) for curr_seq in sequences]
    if rc_truth and return_both:
        rev_dna_seq = [curr_seq.reverse_complement() for curr_seq in dna_seq]
        return rev_dna_seq + dna_seq
    elif rc_truth and not return_both:
        rev_dna_seq = [curr_seq.reverse_complement() for curr_seq in dna_seq]
        return rev_dna_seq
    else:
        return dna_seq


def make_fasta(sequences: List[str], prompts: List[str], ids: List[str], output_file: str) -> None:
    """
    Create a FASTA file from sequences with corresponding IDs and descriptions.

    Args:
        sequences: List of DNA sequences
        prompts: List of sequence descriptions
        ids: List of sequence identifiers
        output_file: Path to save the FASTA file
    Files Generated:
       Creates a FASTA files containing all the generated sequences
       in FASTA format.
    """
    dna_seq_record = [
        SeqRecord(Seq(dna_seq), id=seq_id, description=prompt)
        for dna_seq, seq_id, prompt in zip(sequences, ids, prompts)
    ]
    with open(output_file, "w") as output_handle:
        SeqIO.write(dna_seq_record, output_handle, "fasta")


def sample_model(
    prompt_batches: BatchType,
    model: StripedHyena,
    tokenizer: CharLevelTokenizer,
    file_save_location: str,
    n_tokens: int = 1000,
    temp: float = 0.7,
    top_k: int = 4,
    batched: bool = True,
    n_sample_per_prompt: int = 1,
    device: str = "cuda:0",
    force_prompt_threshold: int = 2,
) -> Tuple[List[str], List[str], List[float], List[str]]:
    """Sample sequences from Evo and save.

    Args:
    prompt_batches: Input sequences organized as batches.
                   For batched=True: List[List[str]] where each inner list contains
                   sequences of the same length.
                   For batched=False: List[str] of individual sequences.

    model: Initialized Evo model for sequence generation.

    tokenizer: CharLevelTokenizer for processing sequences.

    file_save_location: Path where the output CSV will be saved.
                      CSV will contain columns: UUID, Prompt, Generated Sequence, Score

    n_tokens: Maximum number of tokens to generate per sequence.
             Default = 1000.

    temp: Temperature parameter for generation randomness.

    top_k: Number of top tokens to sample from.

    batched: Whether to process prompts in batches (True) or individually (False).

    n_sample_per_prompt: Number of sequences to generate per input prompt.
                       Example: if n_sample_per_prompt=3, each prompt will generate
                       3 different sequences.

    device: Computing device for model execution ('cuda:0', 'cpu', etc).

    force_prompt_threshold: Minimum tokens matching prompt.

    Returns:
         Tuple containing four parallel lists:
         - prompts: Original input sequences
         - sequences: Generated protein sequences
         - scores: Quality scores for generated sequences (as floats)
         - ids: Unique UUIDs assigned to each generation
    Files Generated:
         Creates a CSV file at file_save_location containing:
         - UUID: Unique identifier for each generation
         - Prompt: Original input sequence
         - Generated Sequence: Model output sequence
         - Score: Generation quality score

    """
    print("Starting sample")
    coupledpromptseqscores: List[List[Any]] = []
    extended_header = ["UUID", "Prompt", "Generated Sequence", "Score"]

    if batched:
        for _ in range(n_sample_per_prompt):
            for batch in prompt_batches:
                print("Created batches")
                valid_batch = [seq for seq in batch if isinstance(seq, str) and seq.strip()]
                if not valid_batch:
                    continue
                genseqs, genscores = run_model(
                    valid_batch,
                    model,
                    tokenizer,
                    n_tokens,
                    temp,
                    top_k,
                    batched,
                    device,
                    force_prompt_threshold,
                )
                print("Ran model")

                current_batch = [
                    [uuid.uuid4().hex, prompt, seq, str(score)]
                    for prompt, seq, score in zip(valid_batch, genseqs, genscores)
                    if isinstance(seq, str) and seq.strip() and not math.isnan(float(score))
                ]
                print(current_batch)
                coupledpromptseqscores.extend(current_batch)
    else:
        for i in range(n_sample_per_prompt):
            genseqs, genscores = run_model(
                prompt_batches[i],
                model,
                tokenizer,
                n_tokens,
                temp,
                top_k,
                batched,
                device,
                force_prompt_threshold,
            )
            current_batch = [
                [uuid.uuid4().hex, prompt, seq, str(score)]
                for prompt, seq, score in zip([prompt_batches[i]], genseqs, genscores)
            ]
            coupledpromptseqscores.extend(current_batch)

    filtered_data = []
    for row in coupledpromptseqscores:
        if len(row) < 4:
            continue

        sequence_id, prompt, sequence, score_str = row
        if not isinstance(prompt, str) or not prompt.strip():
            continue
        if not isinstance(sequence, str) or not sequence.strip():
            continue
        try:
            score = float(score_str)
            if math.isnan(score):
                continue
        except (ValueError, TypeError):
            continue
        filtered_data.append([sequence_id, prompt.strip(), sequence.strip(), score])

    # Save filtered results to CSV
    with open(file_save_location, mode="w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(extended_header)
        writer.writerows(filtered_data)

    return read_evo_seqs(filtered_data)


def run_prodigal(
    input_file: str,
    output_file: str,
    output_orf_file: str,
    prodigal_path: str = "/old_home/cirrascale/miniconda3/envs/evo-design/bin/prodigal",
) -> None:
    """
    Run Prodigal for gene prediction on metagenomic sequences.

    Args:
        input_file: Path to input FASTA file containing genomic sequences
        output_file: Path to save predicted protein sequences (amino acids)
        output_orf_file: Path to save predicted coding sequences (nucleotides)
        prodigal_path: Path to Prodigal executable

    Files Generated:
       Creates two FASTA files at output_file (protein sequences) and output_orf_file
       containing the prodigal predicted ORFs and protein sequences called from the
       input sequnces.
    """
    if not Path(prodigal_path).exists():
        raise FileNotFoundError(f"Prodigal not found at: {prodigal_path}")

    cmd = [
        prodigal_path,
        "-i",
        input_file,
        "-a",
        output_file,
        "-d",
        output_orf_file,
        "-p",
        "meta",  # Metagenomics mode
    ]

    subprocess.run(cmd, check=True)


def filter_protein_fasta(
    input_fasta: str,
    output_fasta: str,
    segmasker_path: str,
    min_length: int = 40,
    max_length: int = 1200,
    filter_partial_bool: bool = True,
    segmasker_threshold: float = 0.2,
) -> int:
    """
    Filter protein sequences based on quality criteria. This function applies several filters to identify and remove low-quality or problematic
    protein sequences. The filtering criteria include length constraints, repetitive sequence
    detection, amino acid diversity requirements, and low-complexity region identification
    using NCBI's segmasker tool.

    Args:
        input_fasta: Path to input FASTA file containing protein sequences to filter.
        output_fasta: Path where filtered sequences will be saved.
        segmasker_path: Path to NCBI's segmasker executable.
        min_length: Minimum legnth of protein sequence to be kept.
        max_length: Maximum length of protein sequence to be kept.
        filter_partial_bool: If partial proteins should be filtered out of the generated sequences.
    Returns:
        Number of sequences that passed filters
    Files Generated:
       Creates a FASTA file at output_fasta containing only sequences that passed
       all filters. Original sequence IDs and descriptions are preserved.
    """

    def check_segmasker_installation() -> None:
        """Verify segmasker is installed and provide installation instructions if needed."""
        try:
            result = subprocess.run([segmasker_path, "-version"], capture_output=True, text=True)
            if result.returncode != 0:
                raise FileNotFoundError
        except (FileNotFoundError, subprocess.SubprocessError):
            print("\nError: NCBI segmasker not found or not working properly.")
            print("\nTo install NCBI BLAST+ suite (including segmasker):")
            print("  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
            print("\nAfter installation, provide the correct path to segmasker_path argument")
            raise SystemExit(1)

    def is_segmasked_greater_than_threshold(seq: str, segmasker_threshold: float = 0.1) -> bool:
        """
        Check if sequence has too many low-complexity regions using segmasker.

        Args:
            seq: Protein sequence to analyze
            threshold: Maximum fraction of sequence that can be low-complexity

        Returns:
            True if sequence has too many low-complexity regions
        """
        total_length = len(seq)
        tmp_dir = os.getenv("TMPDIR") or "/dev/shm" if os.path.exists("/dev/shm") else None

        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", dir=tmp_dir) as tmp_file:
            tmp_file.write(f">temp_seq\n{seq}\n".encode())
            tmp_path = tmp_file.name

        try:
            cmd = [
                segmasker_path,
                "-in",
                tmp_path,
                "-outfmt",
                "fasta",
                "-window",
                "15",
                "-locut",
                "1.8",
                "-hicut",
                "3.4",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            return result.stdout.count("X") / total_length > segmasker_threshold
        finally:
            try:
                os.remove(tmp_path)
            except Exception:
                print("segmask failed")
                pass

    def is_highly_repetitive(seq: str, min_repeat_length: int = 3, threshold: float = 0.3) -> bool:
        """
        Check if sequence contains too many repeated elements.

        Args:
            seq: Protein sequence to analyze
            min_repeat_length: Minimum length of repeats to consider
            threshold: Maximum fraction of sequence that can be repetitive

        Returns:
            True if sequence is too repetitive
        """
        seq_len = len(seq)
        seq_array = np.array(list(seq))

        for k in range(min_repeat_length, min_repeat_length + 7):
            kmers = np.lib.stride_tricks.sliding_window_view(seq_array, k)
            kmer_strings = ["".join(kmer) for kmer in kmers]
            if kmer_strings and max(Counter(kmer_strings).values()) * k > seq_len * threshold:
                return True
        return False

    def is_underrepresented_amino_acids(seq: str) -> bool:
        """
        Check if sequence has high proportion of amino acids that occur at very low frequencies.

        Args:
            seq: Protein sequence to analyze

        Returns:
            True if 30% of amino acids are have frequencies under 2, indicative of low sequence complexity
        """
        aa_counts = Counter(seq)
        total_unique = len(aa_counts)
        sorted_counts = sorted(aa_counts.values(), reverse=True)
        num_bottom = max(1, int(0.3 * total_unique))
        return all(count < 2 for count in sorted_counts[-num_bottom:])

    def passes_quality_filters(seq: str, segmasker_threshold: float) -> bool:
        """Apply all filters to a sequence."""
        return not any(
            [
                is_segmasked_greater_than_threshold(seq, segmasker_threshold),
                is_highly_repetitive(seq),
                len(set(seq)) < 12,
                is_underrepresented_amino_acids(seq),
            ]
        )

    def passes_length_and_partial_filters(
        record: str, min_length: int = 40, max_length: int = 1200, filter_partial_bool: bool = False
    ) -> bool:
        """Apply length and partial sequence filters to a record."""
        seq_len = len(record.seq)
        length_ok = min_length <= seq_len <= max_length
        is_complete = True
        if not filter_partial_bool:
            return length_ok
        else:
            is_complete = "partial=00" in record.description
        return length_ok and is_complete

    check_segmasker_installation()

    filtered_records = [
        record
        for record in SeqIO.parse(input_fasta, "fasta")
        if passes_quality_filters(str(record.seq), segmasker_threshold)
        and passes_length_and_partial_filters(record, min_length, max_length, filter_partial_bool)
    ]

    if filtered_records:
        SeqIO.write(filtered_records, output_fasta, "fasta")
    return len(filtered_records)


def run_hmmsearch(
    input_fasta: str, hmm_folder: str, output_csv: str, n_threads: int = 1, incE: float = 0.2
) -> pd.DataFrame:
    """
    Runs all sequences in a FASTA file against all HMMs in a specified folder and records hits in a pandas DataFrame.

    Args:
        input_fasta (str): Path to the input FASTA file with sequences.
        hmm_folder (str): Path to the folder containing HMM files.
        output_csv (str): Path to the output CSV file to save results.
        n_threads (int): Number of parallel threads to use for hmmsearch. Default is 1.

    Returns:
        pd.DataFrame: DataFrame containing the results with columns for sequence, sequence description, HMM, and e-value.
    """
    results = []

    sequences = {record.id: str(record.seq) for record in SeqIO.parse(input_fasta, "fasta")}

    for hmm_file in os.listdir(hmm_folder):
        if hmm_file.endswith(".hmm3"):
            hmm_path = os.path.join(hmm_folder, hmm_file)

            with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".txt") as temp_output:
                temp_output_path = temp_output.name

            try:
                # Run hmmsearch
                subprocess.run(
                    [
                        "hmmsearch",
                        "--domtblout",
                        temp_output_path,  # Output results to a table
                        "--cpu",
                        str(n_threads),  # Number of parallel threads
                        "--incE",
                        str(incE),  # E-value threshold for inclusion
                        hmm_path,  # Input HMM profile
                        input_fasta,  # FASTA with sequences to query against profile
                    ],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )

                # Parse the output and store results
                with open(temp_output_path) as f:
                    for line in f:
                        if not line.startswith("#"):
                            fields = line.split()
                            if len(fields) >= 20:
                                seq_id = fields[0]
                                e_value = float(fields[6])
                                seq_desc = fields[3]
                                seq = sequences.get(seq_id, "")
                                results.append(
                                    {
                                        "Evo Sequence ID": seq_id,
                                        "Sequence Description": seq_desc,
                                        "Amino Acid Sequence": seq,
                                        "HMM": hmm_file,
                                        "E-value": e_value,
                                    }
                                )

            except subprocess.CalledProcessError as e:
                print(f"hmmsearch failed with exit code {e.returncode}")
                print("Error output:", e.stderr.decode())
            except FileNotFoundError:
                print(
                    'hmmsearch is not installed. Please install it using: "conda install -c bioconda hmmer"'
                )

            try:
                os.remove(temp_output_path)
            except OSError as e:
                print(f"Error removing temporary file {temp_output_path}: {e}")

    df_results = pd.DataFrame(results)

    df_results.to_csv(output_csv, index=False)

    return df_results


def get_pfam_hits(
    input_fasta: str, pfam_db_path: str, output_csv: str, n_threads: int = 1, verbose: bool = False
) -> pd.DataFrame:
    """
    Uses HMMER's hmmscan to identify Pfam domain hits in protein sequences.

    Args:
        input_fasta: Path to input FASTA file containing protein sequences
        pfam_db_path: Path to Pfam HMM database file
        output_csv: Path where results will be saved as CSV
        n_threads: Number of CPU threads to use for hmmscan (default: 1)
        verbose: Print additional progress information (default: False)

    Returns:
        pandas.DataFrame containing Pfam hits with columns:
        - target name: Name of the matching Pfam domain
        - accession: Pfam accession number
        - query name: Input sequence identifier
        - amino acid sequence: Full protein sequence
        - E-value: Statistical significance of the match
        - score: Bit score of the match
        - bias: Bias correction term
        - Additional alignment statistics and domain descriptions

    Files Generated:
        - A CSV file at output_csv containing all Pfam hits with detailed
          statistics and sequence information
        - Temporary files (automatically cleaned up):
            * .faa file with protein sequences
            * .txt file with hmmscan output
    """
    aa_sequences = list(SeqIO.parse(input_fasta, "fasta"))
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".faa") as tmp_aa_desc_file:
        records = []
        for aa_seq in aa_sequences:
            new_record = SeqRecord(aa_seq.seq, id=str(aa_seq.seq), description=str(aa_seq.seq))
            records.append(new_record)
        SeqIO.write(records, tmp_aa_desc_file, "fasta")
        aa_desc_file = tmp_aa_desc_file.name

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".txt") as tmp_hmmscan_output:
        hmmscan_output_file = tmp_hmmscan_output.name
        cmd = f"hmmscan --domtblout {hmmscan_output_file} --cpu {n_threads} {pfam_db_path} {aa_desc_file}"
        try:
            subprocess.run(cmd, shell=True)
            if verbose:
                print(f"Running hmmscan with command: {cmd}")
                print(f"Hmmscan output file: {hmmscan_output_file}")
        except subprocess.CalledProcessError as e:
            print("hmmsearch failed with exit code:", e.returncode)
            print("Error output:", e.stderr.decode())
            raise
        except FileNotFoundError:
            sys.stderr.write(
                "The HMM for Pfam is not installed. Please download and install it from: https://www.ebi.ac.uk/interpro/download/pfam/"
            )
            raise

        columns = [
            "target name",
            "accession",
            "query name",
            "amino acid sequence",
            "E-value",
            "score",
            "bias",
            "c-Evalue",
            "i-Evalue",
            "score",
            "bias",
            "hmmfrom",
            "hmm to",
            "alifrom",
            "ali to",
            "envfrom",
            "env to",
            "acc",
            "description of target",
        ]
        pfam = []
        with open(hmmscan_output_file) as f:
            for line in f:
                if not line.startswith("#"):
                    fields = line.split()
                    if len(fields) >= len(columns) - 1:
                        pfam.append(fields[: len(columns) - 1] + [" ".join(fields[len(columns) - 1 :])])

        if verbose:
            if not pfam:
                print("No data parsed from hmmscan output.")
            else:
                print(f"Parsed {len(pfam)} lines from hmmscan output.")

        pfam_hits = pd.DataFrame(pfam, columns=columns)

    seq_descriptions = {record.seq: record.description for record in aa_sequences}
    pfam_hits["sequence description"] = pfam_hits["amino acid sequence"].map(seq_descriptions)

    pfam_hits.to_csv(output_csv, index=False)

    tmp_aa_desc_file.close()
    tmp_hmmscan_output.close()
    for tmp_file in [tmp_aa_desc_file.name, tmp_hmmscan_output.name]:
        try:
            os.remove(tmp_file)
        except OSError as e:
            print(f"Error removing temporary file {tmp_file}: {e}")

    return pfam_hits


def fold_proteins(input_file: str, output_file: str, device: str = "cuda:0") -> pd.DataFrame:
    """
    Predict 3D protein structures using the ESMFold model.

     Args:
        input_file: Path to input FASTA file containing protein sequences
        output_file: Path to save results CSV file
        device: PyTorch device to run predictions on (default: 'cuda:0')

    Returns:
        pandas.DataFrame with columns:
        - Amino Acid Sequence: Input protein sequence
        - Evo Sequence ID: Sequence identifier from FASTA
        - PDB Output: Generated PDB structure string
        - Average pLDDT: Per-residue confidence score (0-100)
        - pTM: Predicted TM-score, indicating global fold quality

    Files Generated:
        CSV file at output_file containing:
        - All columns from the returned DataFrame
        - Can be used to save PDB files later
        - Confidence metrics for structure quality assessment
    """

    esmfold = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
    esmfold = esmfold.to(device)
    esmfold.esm = esmfold.esm.half()
    esmfold_tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    folds = []
    for record in SeqIO.parse(input_file, "fasta"):
        protein_seq = str(record.seq).rstrip("*")  # remove stop codon
        seq_id = record.description
        with torch.inference_mode():
            esmfold_in = esmfold_tokenizer([protein_seq], return_tensors="pt", add_special_tokens=False)
            esmfold_out = esmfold(**esmfold_in.to(device))
            esmfold_out_pdb = esmfold.output_to_pdb(esmfold_out)[0]
            avg_plddt = (esmfold_out["plddt"] * esmfold_out["atom37_atom_exists"]).sum(
                dim=(1, 2)
            ) / esmfold_out["atom37_atom_exists"].sum(dim=(1, 2)).item()
            avg_plddt = avg_plddt.cpu().numpy().item()
            ptm = esmfold_out["ptm"].item() if "ptm" in esmfold_out else None
        folds.append([protein_seq, seq_id, esmfold_out_pdb, avg_plddt, ptm])

    aa_fold_stats = pd.DataFrame(
        folds, columns=["Amino Acid Sequence", "Evo Sequence ID", "PDB Output", "Average pLDDT", "pTM"]
    )
    aa_fold_stats.to_csv(output_file, index=False)
    return aa_fold_stats


def filter_proteins_by_threshold(
    df: pd.DataFrame, output_file: str, plddt_threshold: float = 0.4, ptm_threshold: float = 0.4
) -> pd.DataFrame:
    """
    Filter predicted protein structures based on confidence score thresholds.

    Args:
        df: DataFrame containing ESMFold predictions with 'Average pLDDT' and 'pTM' columns
        output_file: Path to save filtered results CSV
        plddt_threshold: Minimum acceptable pLDDT score (0-1 scale)
        ptm_threshold: Minimum acceptable pTM score (0-1 scale)

    Returns:
        DataFrame containing only structures that pass both confidence thresholds

    Files Generated:
        CSV file at output_file containing the filtered DataFrame entries
    """
    filtered_folds = df[(df["Average pLDDT"] >= plddt_threshold) & (df["pTM"] >= ptm_threshold)]
    filtered_folds.to_csv(output_file, index=False)
    return filtered_folds


def run_foldseek(
    filtered_df: pd.DataFrame, output_file: str, db_path: str, sensitivity: float = 7.5, coverage: float = 0.4
) -> pd.DataFrame:
    """
    Search protein structures against PDB database using Foldseek.

    Args:
        filtered_df: DataFrame containing columns:
            - 'PDB Output': PDB format structure string
            - 'Evo Sequence ID': Unique identifier for the structure
            - 'Amino Acid Sequence': Original protein sequence
        output_file: Path to save search results CSV
        db_path: Path to Foldseek-formatted PDB database
        sensitivity: Search sensitivity
        coverage: Minimum alignment coverage threshold (0-1)

    Returns:
        DataFrame with columns:
        - Amino Acid Sequence: Original query sequence
        - Sequence ID: Query structure identifier
        - Query: Query structure name in results
        - Target: Matching PDB structure identifier
        - Alignment TM-score: Structure similarity score (0-1)
        - LDDT: Local distance difference test score
        - Probability: Match confidence score

    Files Generated:
        CSV file at output_file containing all structure matches with:
        - Alignment quality metrics
        - Target PDB identifiers
        - Confidence scores
    """
    with tempfile.TemporaryDirectory() as tmp_path:
        results = []

        for _, row in filtered_df.iterrows():
            pdb_data = row["PDB Output"]
            pdb_name = row["Evo Sequence ID"].split()[0]

            with tempfile.NamedTemporaryFile(delete=False, dir=tmp_path, suffix=".tsv") as tmp_file:
                output_loc = tmp_file.name
            with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb:
                temp_pdb_path = temp_pdb.name
                temp_pdb.write(pdb_data.encode("utf-8"))

            cmd = [
                "foldseek",
                "easy-search",
                temp_pdb_path,
                db_path,
                output_loc,
                tmp_path,
                "-s",
                str(sensitivity),
                "-c",
                str(coverage),
                "--format-output",
                "query,target,alntmscore,lddt,prob",
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                print(f"Error running command: {e.cmd}")
                print(f"Return code: {e.returncode}")
                print(f"Output: {e.output}")
                print(f"Error output: {e.stderr}")

            os.remove(temp_pdb_path)

            column_names = ["query", "target", "alntmscore", "lddt", "prob"]
            tsv_df = pd.read_csv(output_loc, sep="\t", header=None, names=column_names)

            for _, tsv_row in tsv_df.iterrows():
                results.append(
                    [
                        row["Amino Acid Sequence"],
                        pdb_name,
                        tsv_row["query"],
                        tsv_row["target"],
                        tsv_row["alntmscore"],
                        tsv_row["lddt"],
                        tsv_row["prob"],
                    ]
                )

            try:
                os.remove(output_loc)
            except OSError as e:
                print(f"Error: {e.strerror} - {e.filename}")

        foldseek_df = pd.DataFrame(
            results,
            columns=[
                "Amino Acid Sequence",
                "Sequence ID",
                "Query",
                "Target",
                "Alignment TM-score",
                "LDDT",
                "Probability",
            ],
        )
        foldseek_df.to_csv(output_file, index=False)
        return foldseek_df


def filt_foldseek(foldseek_df: pd.DataFrame, output_csv: str, tm_score_threshold: float = 0.4) -> None:
    """
    Filter Foldseek search results based on TM-score.

    Args:
        foldseek_df: DataFrame containing Foldseek search results with columns:
            - 'Amino Acid Sequence': Protein sequence
            - 'Sequence ID': Sequence identifier
            - 'Alignment TM-score': Structure similarity score
        output_csv: Path to save filtered results
        tm_score_threshold: Minimum TM-score to keep (0-1 scale)
                          Higher values indicate better structural similarity

    Files Generated:
        CSV file at output_csv containing:
        - Amino Acid Sequence: Protein sequence
        - Sequence ID: Sequence identifier
    """
    filtered_data: List[Dict[str, str]] = []
    unique_entries = set()

    for _, row in foldseek_df.iterrows():
        if row["Alignment TM-score"] > tm_score_threshold:
            entry = (row["Amino Acid Sequence"], row["Sequence ID"])
            if entry not in unique_entries:
                unique_entries.add(entry)
                filtered_data.append(
                    {"Amino Acid Sequence": row["Amino Acid Sequence"], "Sequence ID": row["Sequence ID"]}
                )

    with open(output_csv, "w", newline="") as csvfile:
        fieldnames = ["Amino Acid Sequence", "Sequence ID"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for data in filtered_data:
            writer.writerow(data)


def run_mmseqs_search(
    fasta_file: str,
    mmseqs_db: str,
    output_csv: str,
    result_dir: str = "mmseqs_results",
    threads: int = 4,
    sensitivity: float = 4.0,
) -> pd.DataFrame:
    """
    Search protein sequences against a database using MMseqs2.

    Args:
        fasta_file: Path to input FASTA file containing query sequences
        mmseqs_db: Path to MMseqs2-formatted sequence database
        output_csv: Path to save search results CSV
        result_dir: Directory for temporary MMseqs2 files (default: "mmseqs_results")
        threads: Number of CPU threads to use
        sensitivity: Search sensitivity, 1-7 (higher = more sensitive but slower)

    Returns:
        DataFrame containing search results with columns:
        - Query: Query sequence identifier
        - Sequence: Original query sequence
        - Target: Matching sequence identifier
        - Fident: Sequence identity (0-100%)
        - Alnlen: Alignment length
        - Mismatch: Number of mismatches
        - Gapopen: Number of gap openings
        - Qstart/Qend: Query alignment coordinates
        - Tstart/Tend: Target alignment coordinates
        - E-value: Statistical significance
        - Bits: Bit score
        - Qaln/Taln: Aligned sequences

    Files Generated:
        - CSV file at output_csv with all search results
        - Temporary files in result_dir (removed after search)
        - Log file in result_dir/mmseqs_search.log
    """
    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    if not os.path.isdir(mmseqs_db) and not os.path.isfile(mmseqs_db):
        raise FileNotFoundError(f"MMseqs database not found: {mmseqs_db}")

    os.makedirs(result_dir, exist_ok=True)
    mmseqs_out = os.path.join(result_dir, "mmseqs_result.m8")
    log_file = os.path.join(result_dir, "mmseqs_search.log")

    cmd = [
        "mmseqs",
        "easy-search",
        fasta_file,
        mmseqs_db,
        mmseqs_out,
        result_dir,
        "--threads",
        str(threads),
        "-s",
        str(sensitivity),
        "--remove-tmp-files",
        "1",
        "--format-output",
        "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln",
    ]

    print(f"Running MMseqs2 search...")
    try:
        with open(log_file, "w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log, text=True)
    except subprocess.CalledProcessError as e:
        print(f"MMseqs2 search failed with error: {e}")
        raise

    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

    hits = []
    with open(mmseqs_out, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if fields[0] in sequences:
                fields.insert(1, sequences[fields[0]])
                hits.append(fields)

    columns = [
        "Query",
        "Sequence",
        "Target",
        "Fident",
        "Alnlen",
        "Mismatch",
        "Gapopen",
        "Qstart",
        "Qend",
        "Tstart",
        "Tend",
        "E-value",
        "Bits",
        "Qaln",
        "Taln",
    ]

    df = pd.DataFrame(hits, columns=columns)

    numeric_cols = ["Fident", "Alnlen", "Mismatch", "Gapopen", "E-value"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col])

    # Save results
    df.to_csv(output_csv, index=False)

    return df


def align_sequences_mafft(
    input_fasta: str,
    reference_fasta: str,
    output_csv: str,
    output_fasta: str,
    mafft_path: str = "mafft",
    identity_threshold: float = 0.25,
) -> None:
    """
    Align sequences using MAFFT and save results above threshold.

    Args:
        input_fasta: Path to FASTA file containing query sequences
        reference_fasta: Path to FASTA file containing reference sequences
        output_csv: Path to save alignment results CSV
        output_fasta: Path to save filtered sequences FASTA
        mafft_path: Path to MAFFT executable (default: "mafft")
        identity_threshold: Minimum sequence identity to keep (0-1 scale)

    Files Generated:
        - CSV file at output_csv containing columns:
            - Input Sequence Description
            - Input Sequence
            - Best Matching Reference Description
            - Best Matching Reference Sequence
            - Percent Alignment
            - Aligned Input Sequence
            - Aligned Reference Sequence
        - FASTA file at output_fasta containing sequences above threshold
          with original descriptions
    """

    def align_pair(query_record: SeqRecord, ref_record: SeqRecord) -> Tuple[str, str, float]:
        """Align a pair of sequences using MAFFT."""
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

    for filepath in [output_csv, output_fasta]:
        os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)

    query_sequences = list(SeqIO.parse(input_fasta, "fasta"))

    alignment_results = []
    filtered_sequences = []

    for query_record in query_sequences:
        print(f"Processing sequence: {query_record.description}")

        best_match = None
        best_identity = 0
        best_alignment = None

        for ref_record in SeqIO.parse(reference_fasta, "fasta"):
            aligned_query, aligned_ref, identity = align_pair(query_record, ref_record)

            if identity > best_identity:
                best_identity = identity
                best_match = ref_record
                best_alignment = (aligned_query, aligned_ref)

        if best_identity >= identity_threshold:
            alignment_results.append(
                {
                    "Input Sequence Description": query_record.description,
                    "Input Sequence": str(query_record.seq),
                    "Best Matching Reference Description": best_match.description,
                    "Best Matching Reference Sequence": str(best_match.seq),
                    "Percent Alignment": f"{best_identity * 100:.2f}%",
                    "Aligned Input Sequence": best_alignment[0],
                    "Aligned Reference Sequence": best_alignment[1],
                }
            )

            filtered_sequences.append(query_record)

    if alignment_results:
        pd.DataFrame(alignment_results).to_csv(output_csv, index=False)

    if filtered_sequences:
        SeqIO.write(filtered_sequences, output_fasta, "fasta")

    print(f"Processed {len(query_sequences)} sequences")
    print(f"Found {len(filtered_sequences)} sequences above {identity_threshold*100}% identity threshold")
    print(f"Results saved to {output_csv} and {output_fasta}")
