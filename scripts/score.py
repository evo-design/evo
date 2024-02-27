"""
Usage: python -m scripts.score \
           --input-fasta examples/example_seqs.fasta \
           --output-tsv scores.tsv \
           --model-name evo-1-131k-base \
           --device cuda:0

Scores sequences in an input FASTA file according to the joint log-likelihood of the
sequences over all tokens. Outputs these log-likelihood scores to a tab-separated
values file.
"""
import argparse
import pandas as pd

from Bio import SeqIO
from tqdm import tqdm
from evo import Evo, score_sequences


def main():

    # Parse command-line arguments.
    parser = argparse.ArgumentParser(description='Generate sequences using the Evo model.')

    parser.add_argument('--input-fasta', required=True, help='Input FASTA file path')
    parser.add_argument('--output-tsv', required=True, help='Output path to save tab-separated values')
    parser.add_argument('--model-name', type=str, default='evo-1-131k-base', help='Evo model name')
    parser.add_argument('--batch-size', type=int, default=32, help='Number of sequences to evaluate at a time')
    parser.add_argument('--device', type=str, default='cuda:0', help='Device for generation')

    args = parser.parse_args()

    # Load model.

    evo_model = Evo(args.model_name)
    model, tokenizer = evo_model.model, evo_model.tokenizer

    model.to(args.device)
    model.eval()

    # Load sequences.

    seqs = [ str(record.seq) for record in SeqIO.parse(args.input_fasta, 'fasta') ]

    # Score sequences.

    print(f'Scoring {len(seqs)} sequences...')
    scores = []
    for i in tqdm(range(0, len(seqs), args.batch_size)):
        batch_seqs = seqs[i:i + args.batch_size]
        batch_scores = score_sequences(
            batch_seqs,
            model,
            tokenizer,
            device=args.device,
        )
        scores.extend(batch_scores)

    # Save sequences to the output file.

    df = pd.DataFrame({ 'seqs': seqs, 'scores': scores })
    df.to_csv(args.output_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main()