# Copyright (c) Together
# This software is distributed under the terms of the Apache License, Version 2.0
# Author: Michael Poli

import torch
import torch.nn as nn
import torch.nn.functional as F
from einops import rearrange
from torch import Tensor

from .utils import grab_first_if_tuple


class RMSNorm(torch.nn.Module):
    def __init__(self, config):
        super(RMSNorm, self).__init__()
        self.eps, self.hidden_size = config.eps, config.hidden_size
        self.scale = torch.nn.Parameter(torch.ones(self.hidden_size))
        self.register_parameter("scale", self.scale)
        self.scale = self.scale.to(config.params_dtype)
        self.use_flash_rmsnorm = config.get("use_flash_rmsnorm", False)

        if self.use_flash_rmsnorm:
            from flash_attn.ops.rms_norm import rms_norm as rmsnorm_func

            self.rmsnorm_func = rmsnorm_func

    def forward(self, x):
        if self.use_flash_rmsnorm:
            return self.rmsnorm_func(x, self.scale, self.eps)
        else:
            y = x / (x.norm(2, dim=-1, keepdim=True) * self.hidden_size ** (-1.0 / 2) + self.eps)
            return self.scale * y


class ParallelGatedMLP(nn.Module):
    def __init__(
        self,
        config,
    ):
        super().__init__()

        multiple_of = config.get("inner_size_multiple_of", 64)
        self.act_type = config.get("mlp_activation", "gelu")
        if self.act_type == "gelu":
            self.act = F.gelu
        elif self.act_type == "silu":
            self.act = F.silu
        else:
            raise NotImplementedError

        self.multiple_of = multiple_of * config.model_parallel_size

        inner_size = int(2 * config.hidden_size * 4 / 3)
        inner_size = self.multiple_of * ((inner_size + self.multiple_of - 1) // self.multiple_of)
        if config.get("inner_mlp_size", None) is not None:
            inner_size = config.inner_mlp_size

        self.l1 = nn.Linear(
            in_features=config.hidden_size,
            out_features=inner_size,
            bias=False,
        )
        self.l2 = nn.Linear(
            in_features=config.hidden_size,
            out_features=inner_size,
            bias=False,
        )
        self.l3 = nn.Linear(
            in_features=inner_size,
            out_features=config.hidden_size,
            bias=False,
        )

    def forward(self, z):
        z1, z2 = self.l1(z), self.l2(z)
        z1, z2 = grab_first_if_tuple(z1), grab_first_if_tuple(z2)
        y = self.l3(self.act(z1) * z2)
        return grab_first_if_tuple(y)


class Embedding(nn.Module):
    _train_dtype = "bf16"

    def __init__(self, config):
        super().__init__()
        self.word_embeddings = nn.Embedding(config.vocab_size, config.hidden_size, padding_idx=0)

    def embed(self, input_ids, position_ids=None, tokentype_ids=None):
        embeddings = self.word_embeddings(input_ids)
        return embeddings

    def unembed(self, u):
        weight = self.word_embeddings.weight
        return torch.matmul(u, weight)


class VocabParallelEmbedding(nn.Embedding):
    "Adapted from https://github.com/Dao-AILab/flash-attention/blob/main/flash_attn/modules/embedding.py"

    def __init__(self, config):
        vocab_size, process_group, padding_idx = (
            config.vocab_size,
            config.get("process_group", None),
            config.get("padding_idx", None),
        )
        self.process_group = process_group
        if process_group is not None:
            world_size = torch.distributed.get_world_size(process_group)
            if vocab_size % world_size != 0:
                raise ValueError(f"vocab_size ({vocab_size}) must be divisible by " f"world_size ({world_size})")
            if world_size > 1 and padding_idx is not None:
                raise RuntimeError("ParallelEmbedding does not support padding_idx")
        else:
            world_size = 1
        super().__init__(
            vocab_size // world_size,
            embedding_dim=config.hidden_size,
            padding_idx=padding_idx,
        )

    def embed(self, input: Tensor) -> Tensor:
        if self.process_group is None:
            return self.forward(input)
        else:
            rank = torch.distributed.get_rank(self.process_group)
            vocab_size = self.num_embeddings
            vocab_start_index, vocab_end_index = (
                rank * vocab_size,
                (rank + 1) * vocab_size,
            )
            # Create a mask of valid vocab ids (1 means it needs to be masked).
            input_ids_mask = (input < vocab_start_index) | (input >= vocab_end_index)
            input = input - vocab_start_index
            input[input_ids_mask] = 0
            embeddings = self.forward(input)
            embeddings[input_ids_mask] = 0.0
            # Reduce to the global process group
            torch.distributed.all_reduce(embeddings, group=self.process_group)
            return embeddings

    def unembed(self, u: Tensor) -> Tensor:
        if self.process_group is None:
            return u @ self.weight.T
        else:
            raise NotImplementedError
