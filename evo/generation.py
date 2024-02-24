import numpy as np
import sys
from typing import List, Tuple

from .scoring import logits_to_logprobs, prepare_batch
from stripedhyena.generation import Generator
from stripedhyena.model import StripedHyena
from stripedhyena.tokenizer import CharLevelTokenizer


def generate(
    prompt_seqs: List[str],
    model: StripedHyena,
    tokenizer: CharLevelTokenizer,
    n_tokens: int = 100,
    temperature: float = 0.,
    top_k: int = 1,
    top_p: float = 1.,
    batched: bool = True,
    prepend_bos: bool = True,
    cached_generation: bool = False,
    verbose: int = 1,
    device: str = 'cuda:0',
    *args, **kwargs,
) -> Tuple[List[str], List[float]]:
    """
    Performs generation from a list of prompts.
    If all prompts are the same length, this can do batched generation.
    Also supports cached generation for efficient sampling.
    """
    model.eval()

    g = Generator(
        model,
        tokenizer,
        top_k=top_k,
        top_p=top_p,
        temperature=temperature,
    )

    uniform_lengths = all(len(s) == len(prompt_seqs[0]) for s in prompt_seqs)

    if batched and uniform_lengths:
        input_ids_list = [
            prepare_batch(
                prompt_seqs,
                tokenizer,
                prepend_bos=prepend_bos,
                device=device,
            )[0]
        ]
    else:
        if verbose:
            if not uniform_lengths:
                sys.stderr.write('Note: Prompts are of different lengths.\n')
            sys.stderr.write('Note: Will not do batched generation.\n')
        input_ids_list = [
            prepare_batch(
                [ prompt_seq ],
                tokenizer,
                prepend_bos=prepend_bos,
                device=device,
            )[0]
            for prompt_seq in prompt_seqs
        ]

    generated_seqs, generated_scores = [], []
    for input_ids in input_ids_list:
        batch_size = input_ids.shape[0]
        
        output_ids, logits = g.generate(
            input_ids=input_ids,
            num_tokens=n_tokens,
            cached_generation=cached_generation,
            device=device,
            print_generation=True,
            verbose=(verbose > 1),
            stop_at_eos=False,
        )
        if verbose > 1:
            print('input_ids.shape', input_ids.shape)
            print('output_ids.shape', output_ids.shape)
            print('logits.shape', logits.shape)

        generated_seqs_batch = list(tokenizer.detokenize_batch(output_ids))
        assert len(generated_seqs_batch) == batch_size
        generated_seqs += generated_seqs_batch

        logprobs = logits_to_logprobs(logits, output_ids, trim_bos=prepend_bos)
        logprobs = logprobs.float().cpu().numpy()

        generated_scores += [ np.mean(logprobs[idx]) for idx in range(batch_size) ]

    assert len(generated_seqs) == len(generated_scores) == len(prompt_seqs)
    if verbose:
        for seq, score, prompt in zip(generated_seqs, generated_scores, prompt_seqs):
            print(f'Prompt: "{prompt}",\tOutput: "{seq}",\tScore: {score}')

    return generated_seqs, generated_scores
