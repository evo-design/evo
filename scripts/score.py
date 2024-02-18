import argparse
from Bio import SeqIO
import pandas as pd

from evo import Evo, score_sequences


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate sequences using the Evo model.')

    parser.add_argument('--input-fasta', required=True, help='Input FASTA file path')
    parser.add_argument('--output-tsv', required=True, help='Output path to save tab-separated values')
    parser.add_argument('--model-name', type=str, default='evo-1_stripedhyena_pretrained_8k', help='Evo model name')
    parser.add_argument('--batch-size', type=int, default=1000, help='Number of sequences to evaluate at a time')
    parser.add_argument('--device', type=str, default='cuda:0', help='Device for generation')

    args = parser.parse_args()

    evo_model = Evo(args.model_name)
    model, tokenizer = evo_model.model, evo_model.tokenizer

    model.to(args.device)
    model.eval()

    # Load sequences.

    seqs = [ str(record.seq) for record in SeqIO.parse(args.input_fasta, 'fasta') ]

    # Score sequences.

    scores = []
    for i in range(0, len(seqs), args.batch_size):
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
