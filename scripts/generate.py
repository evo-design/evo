"""
Usage: python -m scripts.generate \
           --model-name evo-1-131k-base \
           --prompt ACGT \
           --n-samples 10 \
           --n-tokens 100 \
           --temperature 1. \
           --top-k 4 \
           --device cuda:0

Generates a sequence given a prompt. Also enables the user to specify various basic
sampling parameters.
"""
import argparse

from evo import Evo, generate


def main():

    # Parse command-line arguments.
    parser = argparse.ArgumentParser(description='Generate sequences using the Evo model.')

    parser.add_argument('--model-name', type=str, default='evo-1-131k-base', help='Evo model name')
    parser.add_argument('--prompt', type=str, default='ACGT', help='Prompt for generation')
    parser.add_argument('--n-samples', type=int, default=3, help='Number of sequences to sample at once')
    parser.add_argument('--n-tokens', type=int, default=100, help='Number of tokens to generate')
    parser.add_argument('--temperature', type=float, default=1.0, help='Temperature during sampling')
    parser.add_argument('--top-k', type=int, default=4, help='Top K during sampling')
    parser.add_argument('--top-p', type=float, default=1., help='Top P during sampling')
    parser.add_argument('--cached-generation', type=bool, default=True, help='Use KV caching during generation')
    parser.add_argument('--batched', type=bool, default=True, help='Use batched generation')
    parser.add_argument('--prepend-bos', type=bool, default=False, help='Prepend BOS token')
    parser.add_argument('--device', type=str, default='cuda:0', help='Device for generation')
    parser.add_argument('--verbose', type=int, default=1, help='Verbosity level')

    args = parser.parse_args()

    # Load model.

    evo_model = Evo(args.model_name)
    model, tokenizer = evo_model.model, evo_model.tokenizer

    model.to(args.device)
    model.eval()

    # Sample sequences.
    
    print('Generated sequences:')
    output_seqs, output_scores = generate(
        [ args.prompt ] * args.n_samples,
        model,
        tokenizer,
        n_tokens=args.n_tokens,
        temperature=args.temperature,
        top_k=args.top_k,
        top_p=args.top_p,
        cached_generation=args.cached_generation,
        batched=args.batched,
        prepend_bos=args.prepend_bos,
        device=args.device,
        verbose=args.verbose,
    )


if __name__ == '__main__':
    main()