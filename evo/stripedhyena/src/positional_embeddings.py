"""
Armin Thomas, Jan 2023.  Modified by Eric Nguyen.

Wrappers for linearly interpolated rope embeddings to use inside of MHA layers of Flash Attn.

"""

import torch
import copy
from einops import rearrange
from flash_attn.layers.rotary import RotaryEmbedding
from flash_attn.modules.mha import MHA


# simple wrapper for flash-attn RoPE with linear scaling:
class LinearlyScaledRotaryEmbedding(RotaryEmbedding):
    def __init__(
        self,
        dim: int,
        scaling_factor: float=1.,
        base=10000.0,
        interleaved=False,
        scale_base=None,
        pos_idx_in_fp32=True,
        device=None,
    ):
        super().__init__(
            dim=dim,
            base=base,
            interleaved=interleaved,
            scale_base=scale_base,
            pos_idx_in_fp32=pos_idx_in_fp32,
            device=device
        )
        self._linear_scaling_factor = scaling_factor
    # adpated from: https://github.com/Dao-AILab/flash-attention/blob/43ceab630bc6c27712428da5a33fc9cb5c369d91/flash_attn/layers/rotary.py#L368
    def _update_cos_sin_cache(self, seqlen, device=None, dtype=None):
        # Reset the tables if the sequence length has changed,
        # if we're on a new device (possibly due to tracing for instance),
        # or if we're switching from inference mode to training
        if (
            seqlen > self._seq_len_cached
            or self._cos_cached is None
            or self._cos_cached.device != device
            or self._cos_cached.dtype != dtype
            or (self.training and self._cos_cached.is_inference())
        ):
            self._seq_len_cached = seqlen
            # We want fp32 here, not self.inv_freq.dtype, since the model could be loaded in bf16
            # And the output of arange can be quite large, so bf16 would lose a lot of precision.
            # However, for compatibility reason, we add an option to use the dtype of self.inv_freq.
            if self.pos_idx_in_fp32:
                t = torch.arange(seqlen, device=device, dtype=torch.float32)
                # linear scaling:
                t = t / self._linear_scaling_factor
                # We want fp32 here as well since inv_freq will be multiplied with t, and the output
                # will be large. Having it in bf16 will lose a lot of precision and cause the
                # cos & sin output to change significantly.
                # We want to recompute self.inv_freq if it was not loaded in fp32
                if self.inv_freq.dtype != torch.float32:
                    inv_freq = self._compute_inv_freq(device=device)
                else:
                    inv_freq = self.inv_freq
            else:
                t = torch.arange(seqlen, device=device, dtype=self.inv_freq.dtype)
                # linear scaling:
                t = t / self._linear_scaling_factor
                inv_freq = self.inv_freq
            # Don't do einsum, it converts fp32 to fp16 under AMP
            # freqs = torch.einsum("i,j->ij", t, self.inv_freq)
            freqs = torch.outer(t, inv_freq)
            if self.scale is None:
                self._cos_cached = torch.cos(freqs).to(dtype)
                self._sin_cached = torch.sin(freqs).to(dtype)
            else:
                power = (
                    torch.arange(seqlen, dtype=self.scale.dtype, device=self.scale.device)
                    - seqlen // 2
                ) / self.scale_base
                scale = self.scale.to(device=power.device) ** rearrange(power, "s -> s 1")
                # We want the multiplication by scale to happen in fp32
                self._cos_cached = (torch.cos(freqs) * scale).to(dtype)
                self._sin_cached = (torch.sin(freqs) * scale).to(dtype)
                self._cos_k_cached = (torch.cos(freqs) / scale).to(dtype)
                self._sin_k_cached = (torch.sin(freqs) / scale).to(dtype)

# swap out RoPE of existing mha:
def swap_mha_rope(
    mha,
    new_rope: torch.nn.Module=LinearlyScaledRotaryEmbedding,
    kwargs_new_rope: dict=None
):
    # determine mha dtype and device:
    dtype = mha.Wq.weight.dtype if mha.cross_attn else mha.Wqkv.weight.dtype
    device = mha.Wq.weight.device if mha.cross_attn else mha.Wqkv.weight.device
    # determine RoPE settings:
    kwargs_old_rope = dict(
        dim = mha.rotary_emb.dim,
        base = mha.rotary_emb.base,
        interleaved = mha.rotary_emb.interleaved,
        scale_base = mha.rotary_emb.scale_base,
        pos_idx_in_fp32 = mha.rotary_emb.pos_idx_in_fp32,
        device = mha.rotary_emb.inv_freq.device
    )
    # delete old RoPE:
    del mha.rotary_emb
    # create new RoPE:
    kwargs_new_rope = kwargs_new_rope or {'scaling_factor': 1.0}
    scaled_rope = new_rope(
        **kwargs_new_rope,
        **kwargs_old_rope
    ).to(dtype)
    # attach new RoPE to mha:
    mha.rotary_emb = scaled_rope
    # make new sure RoPE is correctly registered:
    assert isinstance(mha.rotary_emb, new_rope)
    return mha