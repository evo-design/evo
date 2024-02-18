# Copyright (c) Together
# This software is distributed under the terms of the Apache License, Version 2.0
# Author: Michael Poli

from dataclasses import dataclass, field
from typing import Optional

from torch import Tensor


# https://github.com/Dao-AILab/flash-attention/blob/main/flash_attn/utils/generation.py
@dataclass
class InferenceParams:
    """Inference parameters that are passed to the main model in order
    to efficienly calculate and store the context during inference."""

    max_seqlen: int
    max_batch_size: int
    seqlen_offset: int = 0
    batch_size_offset: int = 0
    key_value_memory_dict: dict = field(default_factory=dict)
    lengths_per_sample: Optional[Tensor] = None

    def reset(self, max_seqlen, max_batch_size):
        self.max_seqlen = max_seqlen
        self.max_batch_size = max_batch_size
        self.seqlen_offset = 0
        if self.lengths_per_sample is not None:
            self.lengths_per_sample.zero_()


@dataclass
class RecurrentInferenceParams:
    """Inference parameters passed to blocks with recurrent mode."""

    fir_filter_length: int = 3
    state_dim: int = 16
    # seqlen_offset not used
    seqlen_offset: int = 0
    fir_state_dict: dict = field(default_factory=dict)
    state_dict: dict = field(default_factory=dict)

    def reset(self):
        self.fir_filter_length = 3
        self.state_dim = 16
        self.seqlen_offset = 0
