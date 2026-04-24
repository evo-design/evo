# Memory optimizations for fine-tuning Evo-1 on long DNA contexts

This folder adds one self-contained file, `memory_patches.py`, plus an example
integration in `examples/patch_hf_trainer.py`. No upstream `evo` or
`stripedhyena` source was modified.

The goal: fine-tune `evo-1-131k-base` with LoRA (or optional QLoRA) at
**50,000-token context on 40 GB GPUs** without OOM.

## TL;DR

Apply these five calls around your existing `from_pretrained` + `get_peft_model`:

1. `patch_stripedhyena_flash_attn_compat()`   — before model load
2. `patch_parallel_iir_rfft(enable=True)`     — before model load
3. `patch_compute_filter(model, num_chunks=4)` — after model load
4. `cast_model_to_bf16(model)` (or `cast_non_quantized_model_to_bf16` for QLoRA) — after LoRA
5. `apply_activation_ckpt(model)` — after LoRA, replaces any `gradient_checkpointing_enable()`

See `examples/patch_hf_trainer.py` for the full 30-line example.

## Why QLoRA alone didn't save memory last time

In the earlier `QLoRA + DDP` attempt, QLoRA showed no memory reduction at 50k
context. The reason (speculated):

- QLoRA only quantizes `torch.nn.Linear` weights to 4-bit.
- It does **not** touch:
  - StripedHyena's `poles` / `residues` / `t` (fp32 complex buffers)
  - Long-sequence activations (`x1v`, FFT intermediates, the filter `h`)
  - LoRA adapter weights (fp32 by default)
- At 50k context the peak memory is dominated by a single **~12 GiB transient
  tensor** created inside `compute_filter`. Quantising weights from 14 GiB to
  3.5 GiB helps the base weights, but that 12 GiB ceiling is unaffected, so the
  observed peak barely moved.



**Recommendation:** Apply the long-context patches first. 

## What each patch does

| Patch | What it changes | Expected benefit |
|---|---|---|
| `patch_compute_filter(model, num_chunks=4)` | Splits the `(residues * exp(log_poles * t))` tensor in `ParallelHyenaFilter` along the time dim | **~12 GiB → ~3 GiB** peak for that transient |
| `patch_parallel_iir_rfft(enable=True)` | Uses `rfft(x1v)` instead of `fft(x1v)` in `HyenaInferenceEngine.parallel_iir` | Halves the FFT intermediate for `x1v` (saves a few GiB at 50k) |
| `apply_activation_ckpt(model)` | Non-reentrant activation checkpointing on `AttentionBlock` / `ParallelGatedConvBlock` | Drops activation memory to `~sqrt(L)` scaling |
| `cast_model_to_bf16(model)` | Casts every fp param + buffer (incl. `poles`, `residues`) to bf16 | Halves non-quantised weight memory; avoids the fp32-vs-bf16 FSDP flatten error |
| `cast_non_quantized_model_to_bf16(model)` | Same as above, but preserves bnb 4-bit / 8-bit params | Same benefit under QLoRA |
| `quantize_linear_layers_for_qlora(model, quant_type="nf4")` | Swaps every `nn.Linear` with `bnb.nn.Linear4bit` | 

## Small compatibility shims (also in this PR)

| Patch | Why |
|---|---|
| `patch_stripedhyena_flash_attn_compat()` | `stripedhyena 0.2.2` passes `pos_idx_in_fp32` to flash-attn's `RotaryEmbedding`, which is unexpected in flash-attn ≥ 2.7. Filters kwargs through `inspect.signature`. Safe no-op if not needed. |
| `normalize_model_config_for_peft(model)` | PEFT calls `model.config.to_dict()`. For `AutoModelForCausalLM.from_pretrained` this already works. For Evo's own `load_checkpoint`, `config` is a plain `dict` and this call fails. This shim wraps it. |
| `collect_lora_targets(model)` | Helper that lists both `nn.Linear` and `bnb.nn.Linear4bit` names, if you want to pass an explicit list instead of `target_modules="all-linear"`. |

## Apply order matters

```
patch_stripedhyena_flash_attn_compat()      # before model load
patch_parallel_iir_rfft(enable=True)         # before model load

model = AutoModelForCausalLM.from_pretrained(...)

patch_compute_filter(model, num_chunks=4)    # after load
normalize_model_config_for_peft(model)       # before get_peft_model (safe no-op for HF)

# Optional QLoRA:
# quantize_linear_layers_for_qlora(model, quant_type="nf4")

model = get_peft_model(model, lora_config)   # LoRA sees 4-bit layers if QLoRA ran

# Cast non-quantised tensors to bf16 (must run AFTER LoRA injection):
cast_non_quantized_model_to_bf16(model)      # if QLoRA
# cast_model_to_bf16(model)                  # if no QLoRA

apply_activation_ckpt(model)                 # LAST
```

If `apply_activation_ckpt` runs before the `forward` wrappers you may be using
(e.g. to pop `inputs_embeds` kwargs), the checkpoint wrapper will capture the
un-patched `forward`. Apply any `forward` wrappers first, then
`apply_activation_ckpt`.

**Also remove** the original HuggingFace call
`model.base_model.model.gradient_checkpointing_enable()`. It is a no-op on
StripedHyena and is replaced by `apply_activation_ckpt(model)`.

## Recommended setup for DDP (not FSDP)

For the 7-node × 4-GPU DDP run:

- Stay on **DDP**. FSDP adds a "uniform dtype" constraint that makes LoRA +
  bf16 fragile on StripedHyena, with no real payoff when the only trainable
  parameters are LoRA adapters.
- LoRA is trained; base weights are frozen. Gradient all-reduce cost is tiny.
- `ddp_find_unused_parameters=False` is correct.
- Keep `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`.
- Consider bumping `NCCL_TIMEOUT` to `3600` for 28-GPU eval/save barriers.

## Environment we validated on (single GPU)

These patches were tested end-to-end at 50k context on one A40 (48 GB, sm_86)
with:

| Package | Version |
|---|---|
| Python | 3.12 |
| torch | 2.11.0 + cu128 |
| flash-attn | 2.8.3 (prebuilt wheel) |
| stripedhyena | 0.2.2 |
| transformers | 5.6.1 |
| peft | 0.19.1 |
| accelerate | 1.13.0 |
| datasets | 4.8.4 |
| bitsandbytes | 0.49.2 |

Output adapters: `outputs/single_gpu_50k_qlora_rfft_breakdown/` (local only).

## Assumed target environment (multi-node Slurm DDP)

- Python ≥ 3.10
- torch ≥ 2.3 with CUDA matching your Slurm driver
- stripedhyena ≥ 0.2.0
- peft ≥ 0.10 (for `target_modules="all-linear"` covering `Linear4bit`)
- bitsandbytes ≥ 0.43 (only if you enable QLoRA)
- A100 40 GB or 80 GB, bf16-capable
- HuggingFace cache writable from all ranks (set `HF_HOME` once per job)

If flash-attn in your env is older than 2.7, `patch_stripedhyena_flash_attn_compat`
is a no-op and safe to skip.

## How to verify the patches worked

On any rank-0 process:

```python
import torch
torch.cuda.reset_peak_memory_stats()
# ... run one training step ...
print("peak GiB:", torch.cuda.max_memory_allocated() / 1024**3)
```

Tested peak allocation locally: 36GB
