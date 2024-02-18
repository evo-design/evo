import os
import requests
import torch
from typing import List, Tuple
import yaml

from .stripedhyena.src.utils import dotdict
from .stripedhyena.src.model import StripedHyena
from .stripedhyena.src.tokenizer import CharLevelTokenizer


VALID_MODEL_NAMES = [
    'evo-1_stripedhyena_pretrained_8k',
    'evo-1_stripedhyena_pretrained_131k',
]


class Evo:
    def __init__(self, model_name: str, device: str = 'cuda:0'):
        """
        Loads an Evo model checkpoint given a model name.
        If the checkpoint does not exist, automatically downloads the model to
        `~/.cache/torch/hub/checkpoints`.
        """

        if model_name not in VALID_MODEL_NAMES:
            raise ValueError(
                f'Invalid model name {model_name}. Should be one of: '
                f'{", ".join(VALID_MODEL_NAMES)}.'
            )

        # Download checkpoint.

        home_directory = os.path.expanduser('~')
        download_url = f'https://TODO/checkpoints/{model_name}.pt'
        cache_dir = f'{home_directory}/.cache/torch/hub/checkpoints'
        checkpoint_path = f'{cache_dir}/{model_name}.pt'

        if not os.path.exists(checkpoint_path):
            print(f'Downloading {download_url} to {cache_dir}...')

            if not os.path.exists(cache_dir):
                os.makedirs(cache_dir, exist_ok=True)

            response = requests.get(download_url, stream=True)
            if response.status_code == 200:
                with open(checkpoint_path, 'wb') as f:
                    f.write(response.raw.read())
            else:
                raise Exception(f'Failed to download the file. Status code: {response.status_code}')

        # Load correct config file.

        if model_name == 'evo-1_stripedhyena_pretrained_8k':
            config_path = 'evo/stripedhyena/configs/sh_inference_config_7b.yml'
        elif model_name == 'evo-1_stripedhyena_pretrained_131k':
            config_path = 'evo/stripedhyena/configs/sh_inference_config_7b_rotary_scale_16.yml'
        else:
            raise ValueError(f'Invalid model name {model_name}.')

        # Load model.

        self.model, self.tokenizer = load_checkpoint(
            checkpoint_path,
            model_type='stripedhyena',
            config_path=config_path,
            device=device,
        )
        self.model = self.model.to(device).eval()

        self.device = device


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
