import numpy as np
import sys
import torch
from typing import List, Tuple, Union

from stripedhyena.model import StripedHyena
from stripedhyena.sample import sample
from stripedhyena.tokenizer import CharLevelTokenizer

from .scoring import logits_to_logprobs, prepare_batch


class Generator:
    '''
    Adapted from https://github.com/togethercomputer/stripedhyena.

    Modifications include:
    - `generate()` accepts and returns the recurrent cache state, letting the user
      keep track of it across sampling runs.
    - Able to sample with long token prompts in which the cache is initialized with
      recurrent teacher forcing.
    '''
    def __init__(
        self,
        model: StripedHyena,
        tokenizer: CharLevelTokenizer,
        top_k: int = 50,
        top_p: float = 0.7,
        temperature: float = 1.,
    ):
        self.model = model
        self.tokenizer = tokenizer
        self.top_k = top_k
        self.top_p = top_p
        self.temperature = temperature
        self.untils = ['\n\n']

    def generate(
        self,
        device: str,
        input_string: str = None,
        input_ids: torch.tensor = None,
        num_tokens: int = 32,
        cached_generation: bool = True,
        force_prompt_threshold: int = 128,
        print_generation: bool = True,
        verbose: bool = False,
        skip_special_tokens: bool = False,
        stop_at_eos: bool = True,
        max_seqlen: int = None,
        inference_params_dict: dict = None,
    ) -> Tuple[torch.tensor, torch.tensor, dict]:
        """
        A version of the generate() method that enables passing in and that returns the
        `inference_params_dict` for replaying cached sampling from a given state.
        """
        if isinstance(self.tokenizer.eos, int):
            eos_token_ids = torch.LongTensor([self.tokenizer.eos]).to(device)
        else:
            # is a tensor
            eos_token_ids = self.tokenizer.tokenize(self.tokenizer.eos).to(device)

        if input_ids is None:
            input = self.tokenizer.tokenize(input_string)
            if isinstance(input, list):
                input = torch.LongTensor(input).unsqueeze(0).to(device)
            # is a tensor
            else:
                input = input.unsqueeze(0).to(device)

        else:
            input = input_ids
        x = input

        if max_seqlen is not None:
            x = x[:, -max_seqlen :]

        num_tokens = int(num_tokens)
        batch_size = x.shape[0]

        prompt_length = x.shape[1]
        prompt_forcing = prompt_length > force_prompt_threshold
        if prompt_forcing:
            forced_prompt_length = prompt_length - force_prompt_threshold
            x_force = x[:, force_prompt_threshold:]
            x = x[:, :force_prompt_threshold]
        else:
            forced_prompt_length = 0

        generation = torch.empty(
            x.shape[0],
            num_tokens,
            dtype=torch.long,
            device=x.device,
        )

        scores = torch.empty(
            x.shape[0],
            num_tokens,
            self.tokenizer.vocab_size,
            dtype=torch.float,
            device=x.device,
        )

        if inference_params_dict is not None:
            cached_generation = True
            prefilled = True
            # Ensure that the cached data is loaded on the correct device.
            for key, data in inference_params_dict['mha'].key_value_memory_dict.items():
                inference_params_dict['mha'].key_value_memory_dict[key] = data.to(x.device)
            for key, data in inference_params_dict['hyena'].fir_state_dict.items():
                inference_params_dict['hyena'].fir_state_dict[key] = data.to(x.device)
            for key, data in inference_params_dict['hyena'].state_dict.items():
                inference_params_dict['hyena'].state_dict[key] = data.to(x.device)

        elif cached_generation:
            inference_params_dict = self.model.initialize_inference_params()
            inference_params_dict['mha'].max_batch_size = batch_size
            inference_params_dict['hyena'].max_batch_size = batch_size
            prefilled = False

        if verbose:
            mem_after_tok = torch.cuda.memory_allocated(device=x.device) / 1e9
            print(f'Memory after tokenization: {mem_after_tok} GB')
            print('Starting generation...')
            if input_string is not None:
                print('Prompt: ' + input_string)
            else:
                print(f'Prompt ids: {input_ids} {input_ids.shape}')

        for i in range(forced_prompt_length + num_tokens):
            if prefilled:
                post_prefill = True
            else:
                post_prefill = cached_generation and i > 0

            # prefill then process only the last token
            if post_prefill:
                x = x[:, -1:]
                seqlen_offset = inference_params_dict['mha'].seqlen_offset

                if seqlen_offset == 0:
                    seqlen_offset = input.shape[-1]
                    inference_params_dict['hyena'].seqlen_offset = seqlen_offset
                    inference_params_dict['mha'].seqlen_offset = seqlen_offset
                else:
                    inference_params_dict['mha'].seqlen_offset += 1
                    inference_params_dict['hyena'].seqlen_offset += 1

            # do forward pass with no gradient
            with torch.inference_mode():
                logits, inference_params_dict = self.model(
                    x,
                    inference_params_dict=inference_params_dict,
                )

            last_logits = logits[:, -1]

            if prompt_forcing and i < forced_prompt_length:
                new_idx = x_force[:, i]
            else:
                new_idx = sample(
                    last_logits,
                    top_k=self.top_k,
                    top_p=self.top_p,
                    temperature=self.temperature,
                )

            if stop_at_eos and (generation[0, -2:] == eos_token_ids).all():
                print('Stopping generation at EOS')

            if print_generation and verbose and batch_size == 1:
                print(
                    f'{self.tokenizer.detokenize([new_idx.item()])}',
                    end=' ',
                )

            if prompt_forcing:
                if i >= forced_prompt_length:
                    scores[:, i - forced_prompt_length] = last_logits
                    generation[:, i - forced_prompt_length] = new_idx
            else:
                scores[:, i] = last_logits
                generation[:, i] = new_idx

            if post_prefill:
                x = new_idx[:, None]
            else:
                x = torch.cat([x, new_idx[:, None]], dim=-1)

        if verbose:
            y = self.tokenizer.detokenize_batch(generation[:, : i + 1])

            for until in self.untils:
                if until in y:
                    y = y.split(until)[0]
                    break

            print(f'\nInput: {input_string}, Output: {y}')

            mem_end = torch.cuda.memory_allocated(device=x.device) / 1e9
            print(f'Memory after generation: {mem_end} GB')

        return generation[:, : i + 1], scores[:, : i + 1], inference_params_dict


def generate(
        prompt_seqs: List[str],
        model: StripedHyena,
        tokenizer: CharLevelTokenizer,
        n_tokens: int = 100,
        temperature: float = 0.,
        top_k: int = 1,
        top_p: float = 1.,
        batched: bool = True,
        prepend_bos: bool = False,
        cached_generation: bool = False,
        force_prompt_threshold: int = 128,
        verbose: int = 1,
        device: str = 'cuda:0',
        **kwargs,
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
        
        output_ids, logits, _ = g.generate(
            input_ids=input_ids,
            num_tokens=n_tokens,
            cached_generation=cached_generation,
            force_prompt_threshold=force_prompt_threshold,
            device=device,
            print_generation=(verbose > 1),
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

        logprobs = logits_to_logprobs(logits, output_ids)
        logprobs = logprobs.float().cpu().numpy()

        generated_scores += [ np.mean(logprobs[idx]) for idx in range(batch_size) ]

    assert len(generated_seqs) == len(generated_scores) == len(prompt_seqs)
    if verbose:
        for seq, score, prompt in zip(generated_seqs, generated_scores, prompt_seqs):
            print(f'Prompt: "{prompt}",\tOutput: "{seq}",\tScore: {score}')

    return generated_seqs, generated_scores
