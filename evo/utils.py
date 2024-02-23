import numpy as np
import pandas as pd
from typing import Callable


NTs = 'ACGT'

AAs = 'ACDEFGHIKLMNPQRSTVWY'

AA_TO_CODON = {
    '*': ['TAA', 'TAG', 'TGA'],  # Stop.
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],  # Ala.
    'C': ['TGT', 'TGC'],  # Cys.
    'D': ['GAT', 'GAC'],  # Asp.
    'E': ['GAA', 'GAG'],  # Glu.
    'F': ['TTT', 'TTC'],  # Phe.
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],  # Gly.
    'H': ['CAT', 'CAC'],  # His.
    'I': ['ATT', 'ATC', 'ATA'],  # Ile.
    'K': ['AAA', 'AAG'],  # Lys.
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # Leu.
    'M': ['ATG'],  # Met.
    'N': ['AAT', 'AAC'],  # Asn.
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],  # Pro.
    'Q': ['CAA', 'CAG'],  # Gln.
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # Arg.
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # Ser.
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],  # Thr.
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],  # Val.
    'W': ['TGG'],  # Trp.
    'Y': ['TAT', 'TAC'],  # Tyr.
}

CODON_TO_AA = {
    codon: aa
    for aa, codon_list in AA_TO_CODON.items()
    for codon in codon_list
}

AA_3_TO_1 = {
    "Ala": "A", # Alanine
    "Arg": "R", # Arginine
    "Asn": "N", # Asparagine
    "Asp": "D", # Aspartic acid
    "Cys": "C", # Cysteine
    "Gln": "Q", # Glutamine
    "Glu": "E", # Glutamic acid
    "Gly": "G", # Glycine
    "His": "H", # Histidine
    "Ile": "I", # Isoleucine
    "Leu": "L", # Leucine
    "Lys": "K", # Lysine
    "Met": "M", # Methionine
    "Phe": "F", # Phenylalanine
    "Pro": "P", # Proline
    "Ser": "S", # Serine
    "Thr": "T", # Threonine
    "Trp": "W", # Tryptophan
    "Tyr": "Y", # Tyrosine
    "Val": "V"  # Valine
}


def nucleotide_deep_mutational_scan(sequence: str, ignore_wt: bool = True):
    for idx, wt in enumerate(sequence):
        for mt in NTs:
            if ignore_wt and wt == mt:
                continue
            yield (wt, mt, idx)


def parse_blast_output(output_path: str) -> pd.DataFrame:
    """
    Parses standard blast output with `-outfmt 6`.
    """
    # blast default format output fields.
    blast_table_header = [
        'qacc', 'sacc', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
        'qend', 'sstart', 'send', 'evalue',
    ]

    data = []
    with open(output_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.strip() == '':
                continue
            line = line.strip().split()
            data.append(dict(zip(blast_table_header, line)))

    df = pd.DataFrame(data)
    if len(df) == 0:
        return df
    df['evalue'] = df['evalue'].astype(float)

    return df


def parse_erpin_output(output_path: str, name: str) -> pd.DataFrame:
    """
    Parses ERPIN output. For an example, see `eval/data/example_rho_output.txt`.
    """
    # ERPIN format output fields.
    output_fields = [ 'strand', 'index', 'interval', 'score', 'evalue' ]

    data = []
    with open(output_path, 'r') as f:
        for line in f:
            if line.startswith(f'>{name}'):
                meta = dict(zip(output_fields, f.readline().rstrip().split()))
                sequence = f.readline().rstrip()
                start, end = meta['interval'].split('..')
                data.append([
                    f"{name}_{meta['index']}_{meta['strand']}",
                    sequence,
                    int(start),
                    int(end),
                    '+' if meta['strand'] == 'FW' else '-',
                    meta['score'],
                    float(meta['evalue']),
                ])

    return pd.DataFrame(
        data,
        columns=[
            'id',
            'seq',
            'start',
            'end',
            'strand',
            'score',
            'evalue',
        ],
    )


def parse_hmmsearch_output(output_path: str) -> pd.DataFrame:
    """
    Parses standard hmmsearch output.
    """
    # hmmsearch format output fields.
    hmmsearch_table_header = [
        'target', 'target_acc', 'tlen', 'query', 'query_acc', 'qlen',
        'evalue', 'score', 'bias', 'num', 'of', 'cevalue', 'ievalue',
        'dscore', 'dbias', 'hmm_from', 'hmm_to', 'ali_from', 'ali_to',
        'env_from', 'env_to', 'acc', 'desc',
    ]

    data = []
    with open(output_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            data.append(dict(zip(hmmsearch_table_header, line)))

    return pd.DataFrame(data)


def permutation_test(
    score_func: Callable[[np.array, np.array], float],
    x1: np.array,
    x2: np.array,
    n_permutations: int = 100_000,
) -> float:
    """
    Returns a permutation-based P value. Computes the null distribution by
    shuffling the provided data and recomputing the `score_func`.
    """
    if n_permutations < 1:
        raise ValueError('Number of permutations must be positive.')

    x1, x2 = np.array(x1), np.array(x2)

    observed_score = score_func(x1, x2)

    null_distribution = np.array([
        score_func(x1, np.random.permutation(x2))
        for _ in range(n_permutations)
    ])

    return np.mean(null_distribution >= observed_score)
