"""Example: integrate memory_patches into an existing HuggingFace Trainer script
that fine-tunes `evo-1-131k-base` with LoRA over multi-node DDP.

This is a DROP-IN reference, not a standalone runnable script. Copy the marked
sections into your own training script. See `finetuning/MEMORY_NOTES.md` for
context on each call and for the apply order.

Assumes you are loading the model via HuggingFace:

    AutoModelForCausalLM.from_pretrained(
        "togethercomputer/evo-1-131k-base",
        trust_remote_code=True,
        revision="1.1_fix",
    )
"""

import os
import torch
from transformers import AutoConfig, AutoModelForCausalLM
from peft import LoraConfig, TaskType, get_peft_model

from finetuning.memory_patches import (
    apply_activation_ckpt,
    cast_model_to_bf16,
    cast_non_quantized_model_to_bf16,
    normalize_model_config_for_peft,
    patch_compute_filter,
    patch_parallel_iir_rfft,
    patch_stripedhyena_flash_attn_compat,
    quantize_linear_layers_for_qlora,
)

# --------------------------------------------------------------------------- #
# Toggle: set True to enable optional 4-bit QLoRA on top of the long-context  #
# patches. See MEMORY_NOTES.md for when this is useful.                       #
# --------------------------------------------------------------------------- #
USE_QLORA = False


# 1) Library-level patches. Must run BEFORE model construction.
patch_stripedhyena_flash_attn_compat()
patch_parallel_iir_rfft(enable=True)


# 2) Load the model as you would normally.
config = AutoConfig.from_pretrained(
    "togethercomputer/evo-1-131k-base",
    trust_remote_code=True,
    revision="1.1_fix",
)
config.use_cache = False

model = AutoModelForCausalLM.from_pretrained(
    "togethercomputer/evo-1-131k-base",
    config=config,
    trust_remote_code=True,
    revision="1.1_fix",
    torch_dtype=torch.bfloat16,
)

# If your script fixes the unembed weights from the original checkpoint, keep
# that code here.


# 3) Instance-level patches on the loaded model.
#    num_chunks=4 is the tested default; increase to 8 if you still see OOM.
patch_compute_filter(model, num_chunks=4)
normalize_model_config_for_peft(model)


# 4) Optional 4-bit QLoRA. Must run BEFORE get_peft_model.
if USE_QLORA:
    quantize_linear_layers_for_qlora(
        model,
        quant_type="nf4",
        compute_dtype="bfloat16",
    )


# 5) LoRA injection (unchanged from your original script).
lora_config = LoraConfig(
    r=16,
    lora_alpha=32,
    target_modules="all-linear",
    lora_dropout=0.05,
    task_type=TaskType.CAUSAL_LM,
)
model = get_peft_model(model, lora_config)


# 6) Cast the non-quantised tensors (poles, residues, LoRA, layer norms) to
#    bf16 so every tensor flowing through the forward shares a dtype.
if USE_QLORA:
    cast_non_quantized_model_to_bf16(model)
else:
    cast_model_to_bf16(model)


# 7) (Optional) keep any forward-wrapping you already do, e.g. popping
# `inputs_embeds` kwargs. Place it BEFORE `apply_activation_ckpt` below so the
# checkpoint wrapper captures the patched forward.


# 8) Activation checkpointing on StripedHyena blocks.
#    IMPORTANT: delete the original line `gradient_checkpointing_enable()`
#    from your script. It is a no-op on StripedHyena and is replaced by this.
apply_activation_ckpt(model)


# -- From here, the Trainer / TrainingArguments / trainer.train() are unchanged. --

# Multi-node DDP recommendations (see MEMORY_NOTES.md -> "Recommended setup"):
#   - Do NOT enable FSDP or DeepSpeed.
#   - Keep `ddp_find_unused_parameters=False`.
#   - Export `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`.
#   - Export `NCCL_TIMEOUT=3600` (28-GPU eval barriers can be slow).

if __name__ == "__main__":
    pytorch_alloc = os.environ.get("PYTORCH_CUDA_ALLOC_CONF")
    if "expandable_segments:True" not in (pytorch_alloc or ""):
        print(
            "Hint: export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True "
            "before launching torchrun for smoother long-context allocation."
        )
