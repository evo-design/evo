from contextlib import contextmanager
import io
from functools import partial
import numpy as np
import os
import requests
import sys
import tarfile
import torch
from typing import List, Tuple, Union
import yaml

from .stripedhyena.src.utils import dotdict
from .stripedhyena.src.model import StripedHyena
from .stripedhyena.src.tokenizer import CharLevelTokenizer


class NucleotideModel:
    def __init__(self):
        pass

    def score_batch(self, seqs: List[str], batch_size: int = 1) -> List[float]:
        raise NotImplementedError()


class EvoModel(NucleotideModel):
    def __init__(self, ckpt_path: str, rotary_scale: int = 1, device: str = 'cuda:0'):
        if rotary_scale == 1:
            config_path = 'evo/stripedhyena/configs/sh_inference_config_7b.yml'
        elif rotary_scale == 16:
            config_path = 'evo/stripedhyena/configs/sh_inference_config_7b_rotary_scale_16.yml'
        else:
            raise ValueError(f'Rotary scale {rotary_scale} not supported.')

        self.model, self.tokenizer = load_checkpoint(
            ckpt_path,
            model_type='stripedhyena',
            config_path=config_path,
            device=device,
        )
        self.model = self.model.to(device).eval()

        self.device = device

    def score_sequences(self, seqs: List[str], batch_size: int = 1) -> List[float]:
        from .scoring import score_sequences

        scores = []
        for i in range(0, len(seqs), batch_size):
            batch_seqs = seqs[i:i + batch_size]
            batch_scores = score_sequences(
                batch_seqs,
                self.model,
                self.tokenizer,
                device=self.device,
            )
            scores.extend(batch_scores)

        return scores


def load_model(
        model_name: str,
        device: str = 'cuda:0',
) -> NucleotideModel:
    """
    Loads different Evo checkpoints given the model name.
    """
    if model_name == 'evo-1_stripedhyena_pretrained_8k':
        # TODO: Handle checkpoint better.
        return EvoModel(
            '/scratch/brianhie/dna-gen/checkpoint/evo-1_stripedhyena_pretrained_8k.pt',
            rotary_scale=1,
            device=device,
        )

    elif model_name == 'evo-1_stripedhyena_pretrained_131k':
        # TODO: Handle checkpoint better.
        return EvoModel(
            '/scratch/brianhie/dna-gen/checkpoint/evo-1_stripedhyena_pretrained_131k.pt',
            rotary_scale=16,
            device=device,
        )

    else:
        raise ValueError(f'Invalid model name {model_name}.')


def load_checkpoint(
        ckpt_path: str,
        config_path: str = './evo/stripedhyena/configs/sh_inference_config_7b.yml',
        verbose: int = 0,
        device: str = 'cuda:0',
        **kwargs: dict,
) -> Tuple[StripedHyena, CharLevelTokenizer]:
    """
    Loads a checkpoint from a path and corresponding config.
    """
    global_config = dotdict(yaml.load(open(config_path), Loader=yaml.FullLoader))

    model = StripedHyena(global_config)
    tokenizer = CharLevelTokenizer(512)

    model.load_state_dict(torch.load(ckpt_path), strict=True)

    model.to_bfloat16_except_poles_residues()

    return model, tokenizer
