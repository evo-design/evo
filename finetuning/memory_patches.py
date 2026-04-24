"""Memory-saving patches for fine-tuning Evo-1 on long DNA contexts.

These helpers let you drop the peak GPU memory of a 50k-token Evo-1 training
step from ~40+ GiB to the low 20s on a 40 GiB card, so LoRA / QLoRA fine-tuning
becomes feasible on A100-40GB / A40 / similar.

All functions are runtime monkey-patches. None of them edit `stripedhyena` or
`evo` source files. They match StripedHyena classes by class name, so they
work whether the model was loaded via:

    from transformers import AutoModelForCausalLM
    AutoModelForCausalLM.from_pretrained(..., trust_remote_code=True)

or via Evo's own `load_checkpoint(...)`.

See `finetuning/MEMORY_NOTES.md` for why each patch matters and the order in
which to apply them.
"""

from __future__ import annotations

import inspect
from functools import partial
from types import MethodType, SimpleNamespace
from typing import Optional, Union

import torch
import torch.nn.functional as F


def _is_rank_zero() -> bool:
    try:
        import torch.distributed as dist

        if dist.is_available() and dist.is_initialized():
            return dist.get_rank() == 0
    except Exception:
        pass
    return True


def _log(message: str) -> None:
    if _is_rank_zero():
        print(f"[memory_patches] {message}", flush=True)


# ---------------------------------------------------------------------------
# 1. Long-context computation patches
# ---------------------------------------------------------------------------


def patch_compute_filter(model: torch.nn.Module, num_chunks: int = 4) -> int:
    """Chunk `ParallelHyenaFilter.compute_filter` along the time dimension.

    `compute_filter` materialises a `(residues * exp(log_poles * t))` tensor whose
    peak size scales linearly with sequence length. At 50k context this is the
    single largest transient (~12 GiB in fp32 complex) and is the main cause of
    the long-context OOM errors tracked in MEMORY_NOTES.md.

    Splitting `t` into `num_chunks` along the last dim drops the peak to
    `~12 GiB / num_chunks` with no accuracy change.

    Args:
        model: live model instance (already loaded).
        num_chunks: how many pieces to split `t` into. 4 is the tested default.

    Returns:
        Number of Hyena filters that were patched.
    """
    if num_chunks <= 1:
        return 0

    def compute_filter_chunked(self, L, device):
        self.update_time(L, device)
        filter_dtype = torch.float32
        residues = torch.view_as_complex(self.residues.to(filter_dtype))
        log_poles = torch.view_as_complex(self.poles.to(filter_dtype)).log()
        h_chunks = []
        for t_chunk in torch.chunk(self.t, num_chunks, dim=-1):
            h_chunks.append((residues * (log_poles * t_chunk).exp()).real.sum(1)[None])
        h = torch.cat(h_chunks, dim=-1)
        return h, filter_dtype, log_poles, residues

    patched = 0
    for module in model.modules():
        if module.__class__.__name__ == "ParallelHyenaFilter":
            module.compute_filter = MethodType(compute_filter_chunked, module)
            patched += 1
    _log(f"patched compute_filter on {patched} Hyena filters (chunks={num_chunks}).")
    return patched


def patch_parallel_iir_rfft(enable: bool = True) -> bool:
    """Replace `fft(x1v)` with `rfft(x1v)` inside `HyenaInferenceEngine.parallel_iir`.

    `x1v` is a real tensor, so its full complex FFT is symmetric. Using `rfft`
    halves the size of the intermediate `X` tensor along the frequency axis,
    saving a few GiB at 50k context. Numerics are identical.

    This is a library-level monkey-patch on `stripedhyena.engine.HyenaInferenceEngine`
    and is idempotent: calling it twice is a no-op. Call this BEFORE loading the
    model (the patched method is picked up at forward time regardless, but calling
    it early keeps the code path consistent).

    Returns:
        True if the patch was applied, False otherwise.
    """
    if not enable:
        return False

    try:
        from stripedhyena.engine import HyenaInferenceEngine
    except ImportError:
        _log("stripedhyena not importable; skipping parallel_iir rfft patch.")
        return False

    if getattr(HyenaInferenceEngine, "_rfft_x1v_patched", False):
        return True

    original_parallel_iir = HyenaInferenceEngine.parallel_iir

    def parallel_iir_with_rfft_x1v(
        self,
        z_pre,
        h,
        D,
        L,
        poles,
        residues,
        t,
        dims,
        layer_idx,
        inference_params=None,
        prefill_style="fft",
        fftconv_fn=None,
        padding_mask=None,
        use_flashfft=False,
        column_split_hyena=False,
        long_fir_threshold=None,
    ):
        fft_size = 2 * L
        hidden_size, num_attention_heads, hidden_size_per_attention_head, _, _ = dims
        if column_split_hyena:
            z = z_pre.reshape(
                z_pre.shape[0],
                num_attention_heads,
                3 * hidden_size_per_attention_head,
                z_pre.shape[2],
            )
            x2, x1, v = (
                z[:, :, :hidden_size_per_attention_head],
                z[:, :, hidden_size_per_attention_head : 2 * hidden_size_per_attention_head],
                z[:, :, 2 * hidden_size_per_attention_head :],
            )
            x2, x1, v = (
                x2.reshape(x2.shape[0], -1, x2.shape[-1]),
                x1.reshape(x1.shape[0], -1, x1.shape[-1]),
                v.reshape(v.shape[0], -1, v.shape[-1]),
            )
        else:
            x2, x1, v = z_pre.split([hidden_size, hidden_size, hidden_size], dim=1)

        x1v = x1 * v

        if inference_params is not None and prefill_style == "recurrence":
            y = self.prefill_via_direct_recurrence(
                inference_params=inference_params,
                x1v=x1v,
                L=L,
                poles=poles,
                residues=residues,
            )
        else:
            if use_flashfft and (L % 2) == 0:
                y = fftconv_fn(
                    x1v.to(dtype=torch.bfloat16).contiguous(),
                    h.to(dtype=torch.float32),
                )
                X_s = None
            elif long_fir_threshold is None:
                H = torch.fft.rfft(h.to(dtype=torch.float32), n=fft_size) / fft_size
                if inference_params is None:
                    X = torch.fft.rfft(x1v.to(dtype=torch.float32), n=fft_size)
                    X_s = None
                else:
                    X_s = torch.fft.fft(x1v.to(dtype=torch.float32), n=fft_size)
                    X = X_s[..., : H.shape[-1]]
                if len(z_pre.shape) > 3:
                    H = H.unsqueeze(1)
                y = torch.fft.irfft(X * H, n=fft_size, norm="forward")[..., :L]
            else:
                h = h[0][:, None]
                h = h[..., :long_fir_threshold]
                y = F.conv1d(
                    x1v,
                    h.to(dtype=x1v.dtype),
                    stride=1,
                    groups=x1v.shape[1],
                    padding=h.shape[-1] - 1,
                )[..., :L]

        y = y.to(dtype=x1v.dtype)
        y = (y + x1v * D.unsqueeze(-1)) * x2

        if inference_params is not None:
            if prefill_style == "fft":
                self.prefill_via_modal_fft(
                    inference_params=inference_params,
                    x1v=x1v,
                    X_s=X_s,
                    L=L,
                    t=t,
                    poles=poles,
                    dims=dims,
                    layer_idx=layer_idx,
                    use_flashfft=use_flashfft,
                    fftconv_fn=fftconv_fn,
                )
            elif prefill_style != "recurrence":
                raise NotImplementedError

        return y.permute(0, 2, 1)

    HyenaInferenceEngine.parallel_iir = parallel_iir_with_rfft_x1v
    HyenaInferenceEngine._rfft_x1v_original_parallel_iir = original_parallel_iir
    HyenaInferenceEngine._rfft_x1v_patched = True
    _log("patched HyenaInferenceEngine.parallel_iir to use rfft(x1v).")
    return True


def apply_activation_ckpt(model: torch.nn.Module) -> int:
    """Wrap StripedHyena blocks (`AttentionBlock`, `ParallelGatedConvBlock`) with
    non-reentrant activation checkpointing.

    Must be NON-REENTRANT (reentrant checkpointing breaks PEFT gradient flow).
    Apply AFTER LoRA has been injected.

    Replaces the standard HF `gradient_checkpointing_enable()` call, which is a
    no-op on StripedHyena because it is not a HuggingFace-native architecture.

    Returns:
        Number of distinct block classes wrapped.
    """
    from torch.distributed.algorithms._checkpoint.checkpoint_wrapper import (
        CheckpointImpl,
        apply_activation_checkpointing,
        checkpoint_wrapper,
    )

    block_types = tuple(
        {
            type(module)
            for module in model.modules()
            if module.__class__.__name__ in {"AttentionBlock", "ParallelGatedConvBlock"}
        }
    )
    if not block_types:
        _log("no StripedHyena blocks found; skipping activation checkpointing.")
        return 0

    wrapper = partial(checkpoint_wrapper, checkpoint_impl=CheckpointImpl.NO_REENTRANT)
    apply_activation_checkpointing(
        model,
        checkpoint_wrapper_fn=wrapper,
        check_fn=lambda submodule: isinstance(submodule, block_types),
    )
    _log(f"activation checkpointing wrapped {len(block_types)} block type(s).")
    return len(block_types)


# ---------------------------------------------------------------------------
# 2. Compatibility shims (small diffs that are needed for training to start)
# ---------------------------------------------------------------------------


def normalize_model_config_for_peft(model: torch.nn.Module) -> None:
    """Ensure `model.config` exposes a `.to_dict()` method so PEFT can read it.

    StripedHyena models loaded via Evo's `load_checkpoint` keep `model.config`
    as a plain Python `dict`, but `peft.get_peft_model` calls `config.to_dict()`.
    Models loaded via HuggingFace `AutoModelForCausalLM.from_pretrained` already
    have a proper `PretrainedConfig`, so this function is a safe no-op for them.

    Call BEFORE `get_peft_model(model, lora_config)`.
    """
    config = getattr(model, "config", None)
    to_dict = getattr(config, "to_dict", None)
    if callable(to_dict):
        return

    if isinstance(config, dict):
        config_dict = dict(config)
    elif hasattr(config, "items"):
        config_dict = dict(config.items())
    elif config is None:
        config_dict = {"model_type": "stripedhyena"}
    else:
        config_dict = vars(config).copy()

    class _PeftCompatibleConfig(SimpleNamespace):
        def to_dict(self):
            return vars(self).copy()

    model.config = _PeftCompatibleConfig(**config_dict)
    _log("wrapped model.config with a PEFT-compatible to_dict() view.")


def patch_stripedhyena_flash_attn_compat() -> bool:
    """Make stripedhyena's rotary-embedding swap robust to newer flash-attn.

    In `stripedhyena==0.2.2`, `swap_mha_rope` passes kwargs to flash-attn's
    `RotaryEmbedding.__init__` that no longer exist in newer flash-attn releases
    (>= 2.7). That produces `TypeError: __init__() got an unexpected keyword
    argument 'pos_idx_in_fp32'` at model construction time.

    This function replaces `swap_mha_rope` with a version that filters kwargs
    through `inspect.signature`, so construction works on both old and new
    flash-attn.

    Only needed if you see that TypeError. Safe to call either way.

    Returns:
        True if patched, False if deps missing.
    """
    try:
        from einops import rearrange
        from flash_attn.layers.rotary import RotaryEmbedding
        import stripedhyena.model as sh_model
        import stripedhyena.positional_embeddings as sh_positional_embeddings
    except ImportError:
        _log("flash_attn or stripedhyena not importable; skipping rope compat patch.")
        return False

    class CompatibleLinearlyScaledRotaryEmbedding(RotaryEmbedding):
        def __init__(
            self,
            dim: int,
            scaling_factor: float = 1.0,
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
                device=device,
            )
            self.pos_idx_in_fp32 = pos_idx_in_fp32
            self._linear_scaling_factor = scaling_factor

        def _update_cos_sin_cache(self, seqlen, device=None, dtype=None):
            if (
                seqlen > self._seq_len_cached
                or self._cos_cached is None
                or self._cos_cached.device != device
                or self._cos_cached.dtype != dtype
                or (self.training and self._cos_cached.is_inference())
            ):
                self._seq_len_cached = seqlen
                if self.pos_idx_in_fp32:
                    t = torch.arange(seqlen, device=device, dtype=torch.float32)
                    t = t / self._linear_scaling_factor
                    if self.inv_freq.dtype != torch.float32:
                        inv_freq = self._compute_inv_freq(device=device)
                    else:
                        inv_freq = self.inv_freq
                else:
                    t = torch.arange(seqlen, device=device, dtype=self.inv_freq.dtype)
                    t = t / self._linear_scaling_factor
                    inv_freq = self.inv_freq
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
                    self._cos_cached = (torch.cos(freqs) * scale).to(dtype)
                    self._sin_cached = (torch.sin(freqs) * scale).to(dtype)
                    self._cos_k_cached = (torch.cos(freqs) / scale).to(dtype)
                    self._sin_k_cached = (torch.sin(freqs) / scale).to(dtype)

    def compatible_swap_mha_rope(mha, new_rope=CompatibleLinearlyScaledRotaryEmbedding, kwargs_new_rope=None):
        dtype = mha.Wq.weight.dtype if mha.cross_attn else mha.Wqkv.weight.dtype
        old_rope = mha.rotary_emb
        kwargs_old_rope = {
            "dim": old_rope.dim,
            "base": old_rope.base,
            "interleaved": old_rope.interleaved,
        }
        if hasattr(old_rope, "scale_base"):
            kwargs_old_rope["scale_base"] = old_rope.scale_base
        if hasattr(old_rope, "pos_idx_in_fp32"):
            kwargs_old_rope["pos_idx_in_fp32"] = old_rope.pos_idx_in_fp32
        if hasattr(old_rope, "inv_freq"):
            kwargs_old_rope["device"] = old_rope.inv_freq.device

        kwargs_new_rope = kwargs_new_rope or {"scaling_factor": 1.0}
        allowed = set(inspect.signature(new_rope.__init__).parameters)
        filtered_old = {k: v for k, v in kwargs_old_rope.items() if k in allowed}
        filtered_new = {k: v for k, v in kwargs_new_rope.items() if k in allowed}

        del mha.rotary_emb
        mha.rotary_emb = new_rope(**filtered_new, **filtered_old).to(dtype)
        assert isinstance(mha.rotary_emb, new_rope)
        return mha

    sh_positional_embeddings.LinearlyScaledRotaryEmbedding = CompatibleLinearlyScaledRotaryEmbedding
    sh_positional_embeddings.swap_mha_rope = compatible_swap_mha_rope
    sh_model.swap_mha_rope = compatible_swap_mha_rope
    _log("patched stripedhyena swap_mha_rope for newer flash_attn.")
    return True


# ---------------------------------------------------------------------------
# 3. QLoRA (optional, stackable on top of the long-context patches)
# ---------------------------------------------------------------------------


def _is_quantized_parameter(parameter: torch.nn.Parameter) -> bool:
    return parameter.__class__.__name__ in {"Params4bit", "Int8Params"}


def _is_bnb_linear(module: torch.nn.Module) -> bool:
    try:
        import bitsandbytes as bnb
    except ImportError:
        return False
    return isinstance(module, (bnb.nn.Linear4bit, bnb.nn.Linear8bitLt))


def cast_model_to_bf16(model: torch.nn.Module) -> None:
    """Cast every floating-point parameter and buffer to bfloat16.

    Fixes the FSDP / dtype-mismatch errors in MEMORY_NOTES.md: StripedHyena
    keeps `poles` / `residues` / some buffers in fp32 by default, which
    conflicts with bf16 base weights during forward and during FSDP flatten.

    Use this when you are NOT using QLoRA. If you are using QLoRA, call
    `cast_non_quantized_model_to_bf16` instead.
    """
    for parameter in model.parameters():
        if parameter.is_floating_point():
            parameter.data = parameter.data.to(torch.bfloat16)
    for buffer_name, buffer in model.named_buffers():
        if buffer.is_floating_point():
            module_name = buffer_name.rsplit(".", 1)[0] if "." in buffer_name else ""
            owning = model.get_submodule(module_name) if module_name else model
            owning._buffers[buffer_name.split(".")[-1]] = buffer.to(torch.bfloat16)


def cast_non_quantized_model_to_bf16(model: torch.nn.Module) -> None:
    """Same as `cast_model_to_bf16`, but leaves bitsandbytes 4-bit/8-bit params alone.

    Use this AFTER `quantize_linear_layers_for_qlora` and AFTER LoRA injection,
    so that:
      - 4-bit weights stay 4-bit
      - poles / residues / LoRA A/B / layer norms all become bf16
    """
    for parameter in model.parameters():
        if parameter.is_floating_point() and not _is_quantized_parameter(parameter):
            parameter.data = parameter.data.to(torch.bfloat16)
    for buffer_name, buffer in model.named_buffers():
        if buffer.is_floating_point():
            module_name = buffer_name.rsplit(".", 1)[0] if "." in buffer_name else ""
            owning = model.get_submodule(module_name) if module_name else model
            owning._buffers[buffer_name.split(".")[-1]] = buffer.to(torch.bfloat16)


def quantize_linear_layers_for_qlora(
    model: torch.nn.Module,
    quant_type: str = "nf4",
    compute_dtype: Union[str, torch.dtype] = "bfloat16",
    compress_statistics: bool = True,
) -> int:
    """Replace every `torch.nn.Linear` in `model` with `bnb.nn.Linear4bit`.

    Must be called BEFORE `get_peft_model(...)`. LoRA adapters will then be
    attached on top of the 4-bit weights.

    Args:
        model: live model instance.
        quant_type: "nf4" (recommended) or "fp4".
        compute_dtype: dtype used for dequantised matmul. bfloat16 matches the
            rest of the model; keep it at "bfloat16" unless you know why.
        compress_statistics: passes through to bnb.

    Returns:
        Number of Linear layers that were replaced.
    """
    try:
        import bitsandbytes as bnb
    except ImportError as exc:
        raise RuntimeError(
            "bitsandbytes is required for quantize_linear_layers_for_qlora. "
            "Install it with `pip install bitsandbytes`."
        ) from exc

    if isinstance(compute_dtype, str):
        compute_dtype = {
            "bfloat16": torch.bfloat16,
            "float16": torch.float16,
            "float32": torch.float32,
        }[compute_dtype]

    replaced = 0

    def replace(parent):
        nonlocal replaced
        for child_name, child in list(parent.named_children()):
            if isinstance(child, torch.nn.Linear):
                quantized = bnb.nn.Linear4bit(
                    child.in_features,
                    child.out_features,
                    bias=child.bias is not None,
                    compute_dtype=compute_dtype,
                    quant_type=quant_type,
                    compress_statistics=compress_statistics,
                )
                quantized.load_state_dict(child.state_dict())
                quantized.requires_grad_(False)
                setattr(parent, child_name, quantized)
                replaced += 1
            else:
                replace(child)

    replace(model)
    if replaced == 0:
        raise ValueError("No torch.nn.Linear modules found to quantize for QLoRA.")
    _log(f"quantized {replaced} Linear layers to 4-bit ({quant_type}, compute_dtype={compute_dtype}).")
    return replaced


def collect_lora_targets(model: torch.nn.Module) -> list[str]:
    """Return the dotted module names of every LoRA-compatible layer in `model`.

    Useful when you want to pass an explicit `target_modules=[...]` to
    `LoraConfig` instead of the shorthand `"all-linear"`. Covers both plain
    `torch.nn.Linear` and `bitsandbytes.nn.Linear4bit / Linear8bitLt`.
    """
    targets: list[str] = []
    for name, module in model.named_modules():
        if isinstance(module, torch.nn.Linear) or _is_bnb_linear(module):
            targets.append(name)
    if not targets:
        raise ValueError("No supported Linear modules found for LoRA injection.")
    return sorted(set(targets))


__all__ = [
    "patch_compute_filter",
    "patch_parallel_iir_rfft",
    "apply_activation_ckpt",
    "normalize_model_config_for_peft",
    "patch_stripedhyena_flash_attn_compat",
    "cast_model_to_bf16",
    "cast_non_quantized_model_to_bf16",
    "quantize_linear_layers_for_qlora",
    "collect_lora_targets",
]
