import numpy as np
import torch
from typing import List, Tuple

from stripedhyena.model import StripedHyena
from stripedhyena.tokenizer import CharLevelTokenizer


def prepare_batch(
    seqs: List[str],
    tokenizer: CharLevelTokenizer,
    prepend_bos: bool = True,
    device: str = 'cuda:0'
) -> Tuple[torch.Tensor, List[int]]:
    """
    Takes in a list of sequences, tokenizes them, and puts them in a tensor batch.
    If the sequences have differing lengths, then pad up to the maximum sequence length.
    """
    seq_lengths = [ len(seq) for seq in seqs ]
    max_seq_length = max(seq_lengths)

    input_ids = []
    for seq in seqs:
        padding = [tokenizer.pad_id] * (max_seq_length - len(seq))
        input_ids.append(
            torch.tensor(
                ([tokenizer.eod_id] * int(prepend_bos)) + tokenizer.tokenize(seq) + padding,
                dtype=torch.long,
            ).to(device).unsqueeze(0)
        )
    input_ids = torch.cat(input_ids, dim=0)

    return input_ids, seq_lengths


def logits_to_logprobs(
    logits: torch.Tensor,
    input_ids: torch.Tensor,
    trim_bos: bool = True,
) -> torch.Tensor:
    """
    Takes in a tensor of logits of dimension (batch, length, vocab).
    Computes the log-likelihoods using a softmax along the vocab dimension.
    Uses the `input_ids` to index into the log-likelihoods and returns the likelihood
    of the provided sequence at each position with dimension (batch, length).
    """
    softmax_logprobs = torch.log_softmax(logits, dim=-1)
    if trim_bos:
        softmax_logprobs = softmax_logprobs[:, :-1] # Remove last prediction.
        input_ids = input_ids[:, 1:] # Trim BOS added by tokenizer.
    assert(softmax_logprobs.shape[1] == input_ids.shape[1])

    logprobs = torch.gather(
        softmax_logprobs,       # Gather likelihoods...
        2,                      # along the vocab dimension...
        input_ids.unsqueeze(-1) # using the token ids to index.
    ).squeeze(-1)

    return logprobs


def score_sequences(
    seqs: List[str],
    model: StripedHyena,
    tokenizer: CharLevelTokenizer,
    reduce_method: str = 'mean',
    device: str = 'cuda:0',
) -> List[float]:
    """
    Computes the model log-likelihood scores for sequences in `seqs`.
    Uses `reduce_method` to take the mean or sum across the likelihoods at each 
    position (default: `'mean'`).

    Returns a list of scalar scores corresponding to the reduced log-likelihoods for
    each sequence.    
    """
    input_ids, seq_lengths = prepare_batch(seqs, tokenizer, device=device, prepend_bos=True)
    assert(len(seq_lengths) == input_ids.shape[0])

    with torch.inference_mode():
        logits, _ = model(input_ids) # (batch, length, vocab)

    logprobs = logits_to_logprobs(logits, input_ids, trim_bos=True)
    logprobs = logprobs.float().cpu().numpy()

    if reduce_method == 'mean':
        reduce_func = np.mean
    elif reduce_method == 'sum':
        reduce_func = np.sum
    else:
        raise ValueError(f'Invalid reduce_method {reduce_method}')

    return [
        reduce_func(logprobs[idx][:seq_lengths[idx]])
        for idx in range(len(seq_lengths))
    ]


def positional_entropies(
    seqs: List[str],
    model: StripedHyena,
    tokenizer: CharLevelTokenizer,
    device: str = 'cuda:0',
) -> List[np.array]:
    """
    Computes the positional entropies for sequences in `seqs`.

    Returns a list of arrays, where each array is the same length as the
    corresponding sequence length. Each array contains the per-position entropy
    across the vocab dimension.
    """
    input_ids, seq_lengths = prepare_batch(seqs, tokenizer, device=device, prepend_bos=True)
    assert(len(seq_lengths) == input_ids.shape[0])

    with torch.inference_mode():
        logits, _ = model(input_ids) # (batch, length, vocab)
    
    # Tokenizer prepends BOS, remember to remove last prediction.
    softmax_logprobs = torch.log_softmax(logits, dim=-1)[:, :-1]

    entropies = -torch.sum(torch.exp(softmax_logprobs) * softmax_logprobs, dim=-1)
    entropies = entropies.float().cpu().numpy()

    sequence_entropies = [
        entropies[idx][:seq_lengths[idx]] for idx in range(len(seq_lengths))
    ]
    assert all(
        len(seq) == len(entropy) for seq, entropy in zip(seqs, sequence_entropies)
    )

    return sequence_entropies
