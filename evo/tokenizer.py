"""Numpy 2.0-compatible CharLevelTokenizer.

Drop-in replacement for stripedhyena.tokenizer.CharLevelTokenizer that uses
np.frombuffer instead of the removed np.fromstring (binary mode).
"""
from typing import List, Union

import numpy as np
import torch


class CharLevelTokenizer:
    """Character-level tokenizer that maps text to/from byte values."""

    def __init__(self, vocab_size: int = 512):
        self.name = "CharLevelTokenizer"
        self._vocab_size = vocab_size
        self.eod_id = 0
        self.eos_id = 0
        self.pad_id = 1

    def clamp(self, n):
        return max(32, min(n, self.vocab_size))

    @property
    def vocab_size(self):
        return self._vocab_size

    @property
    def eod(self):
        return self.eod_id

    @property
    def eos(self):
        return self.eod_id

    def decode_token(self, token: int):
        return str(chr(self.clamp(token)))

    def tokenize(self, text: str):
        return list(np.frombuffer(text.encode(), dtype=np.uint8))

    def tokenize_batch(self, text_batch: Union[List[str], str]):
        if isinstance(text_batch, list):
            return [self.tokenize(s) for s in text_batch]
        else:
            return self.tokenize(text_batch)

    def detokenize(self, token_ids):
        return "".join(list(map(self.decode_token, token_ids)))

    def detokenize_batch(self, token_ids: Union[List[str], str]):
        if isinstance(token_ids, list):
            return [self.detokenize(s) for s in token_ids]
        elif isinstance(token_ids, torch.Tensor):
            return [self.detokenize(s) for s in token_ids.tolist()]
        else:
            return self.detokenize(token_ids)
