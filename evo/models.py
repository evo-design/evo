import json
import os
import pkgutil
import re

import yaml

from stripedhyena.utils import dotdict
from stripedhyena.model import StripedHyena
from .tokenizer import CharLevelTokenizer


MODEL_NAMES = [
    'evo-1.5-8k-base',
    'evo-1-8k-base',
    'evo-1-131k-base',
    'evo-1-8k-crispr',
    'evo-1-8k-transposon',
]

class Evo:
    def __init__(self, model_name: str = MODEL_NAMES[1], device: str = None):
        """
        Loads an Evo model checkpoint given a model name.
        If the checkpoint does not exist, we automatically download it from HuggingFace.
        """
        self.device = device

        # Check model name.

        if model_name not in MODEL_NAMES:
            raise ValueError(
                f'Invalid model name {model_name}. Should be one of: '
                f'{", ".join(MODEL_NAMES)}.'
            )

        # Assign config path.

        if model_name == 'evo-1-8k-base' or \
           model_name == 'evo-1-8k-crispr' or \
           model_name == 'evo-1-8k-transposon' or \
           model_name == 'evo-1.5-8k-base':
            config_path = 'configs/evo-1-8k-base_inference.yml'
        elif model_name == 'evo-1-131k-base':
            config_path = 'configs/evo-1-131k-base_inference.yml'
        else:
            raise ValueError(
                f'Invalid model name {model_name}. Should be one of: '
                f'{", ".join(MODEL_NAMES)}.'
            )

        # Load model.

        self.model = load_checkpoint(
            model_name=model_name,
            config_path=config_path,
            device=self.device
        )

        # Load tokenizer.

        self.tokenizer = CharLevelTokenizer(512)

        
HF_MODEL_NAME_MAP = {
    'evo-1.5-8k-base': 'evo-design/evo-1.5-8k-base',
    'evo-1-8k-base': 'togethercomputer/evo-1-8k-base',
    'evo-1-131k-base': 'togethercomputer/evo-1-131k-base',
    'evo-1-8k-crispr': 'LongSafari/evo-1-8k-crispr',
    'evo-1-8k-transposon': 'LongSafari/evo-1-8k-transposon',
}

def load_checkpoint(
    model_name: str = MODEL_NAMES[1],
    config_path: str = 'evo/configs/evo-1-131k-base_inference.yml',
    device: str = None,
    *args, **kwargs
):
    """
    Load checkpoint from HuggingFace and place it into SH model.

    Loads safetensors directly via huggingface_hub + safetensors to avoid
    compatibility issues with different transformers versions (e.g., v5
    breaks tied-weight handling for evo-1 models).
    """
    from huggingface_hub import snapshot_download
    from safetensors.torch import load_file

    # Map model name to HuggingFace model name.

    hf_model_name = HF_MODEL_NAME_MAP[model_name]
    revision = '1.1_fix' if re.match(r'evo-1-.*-base', model_name) else 'main'

    # Download / locate cached model files.

    model_dir = snapshot_download(
        hf_model_name,
        revision=revision,
    )

    # Load safetensors weights (multi-shard or single-shard).

    index_path = os.path.join(model_dir, 'model.safetensors.index.json')
    single_path = os.path.join(model_dir, 'model.safetensors')

    raw_state_dict = {}
    if os.path.exists(index_path):
        with open(index_path) as f:
            index = json.load(f)
        shard_files = sorted(set(index['weight_map'].values()))
        for shard_file in shard_files:
            shard_path = os.path.join(model_dir, shard_file)
            raw_state_dict.update(load_file(shard_path))
    elif os.path.exists(single_path):
        raw_state_dict = load_file(single_path)
    else:
        raise FileNotFoundError(
            f'No safetensors files found in {model_dir}. '
            f'Expected model.safetensors.index.json or model.safetensors.'
        )

    # Strip 'backbone.' prefix from keys (HF checkpoint wraps the backbone).

    state_dict = {}
    for key, value in raw_state_dict.items():
        if key.startswith('backbone.'):
            state_dict[key[len('backbone.'):]] = value
        else:
            state_dict[key] = value
    del raw_state_dict

    # Handle tied embeddings: if unembed.weight is missing, copy from
    # embedding_layer.weight. This is needed for evo-1 models whose
    # checkpoints only store the embedding weight once.

    if 'unembed.weight' not in state_dict and 'embedding_layer.weight' in state_dict:
        state_dict['unembed.weight'] = state_dict['embedding_layer.weight']

    # Load SH config.

    config = yaml.safe_load(pkgutil.get_data(__name__, config_path))
    global_config = dotdict(config, Loader=yaml.FullLoader)

    # Load SH Model.

    model = StripedHyena(global_config)
    model.load_state_dict(state_dict, strict=True)
    model.to_bfloat16_except_poles_residues()
    if device is not None:
        model = model.to(device)

    return model
