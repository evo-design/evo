# Evo: DNA foundation modeling from molecular to genome scale

## Setup

### Requirements

Evo uses (FlashAttention-2)[https://github.com/Dao-AILab/flash-attention], which may not work on all GPU architectures.
Please consult the (FlashAttention GitHub repository)[https://github.com/Dao-AILab/flash-attention#installation-and-features] for the current list of supported GPUs.

Evo also uses PyTorch. Make sure the correct (PyTorch version is installed)[https://pytorch.org/] on your system.

### Installation

You can install Evo using `pip`
```bash
pip install evo-model
```
or directly from the GitHub source
```bash
git clone https://github.com/evo-design/evo.git
cd evo/
pip install .
```

## Usage

You can download Devo and use it locally through the Python API. For example:
```python
from evo import Evo
import torch

device = 'cuda:0'

evo_model = Evo('evo-1_stripedhyena_pretrained_8k')
model, tokenizer = evo_model.model, evo_model.tokenizer
model.to(device)
model.eval()

sequence = 'ACGT'
input_ids = torch.tensor(
    tokenizer.tokenize(sequence),
    dtype=torch.int,
).to(device).unsqueeze(0)
logits, _ = model(input_ids) # (batch, length, vocab)

print('Logits: ', logits)
print('Shape (batch, length, vocab): ', logits.shape)
```
Examples of batched inference can be found in [`scripts/example_inference.py`](scripts/example_inference.py).

We provide an example script for how to prompt the model and sample a set of sequences given the prompt.
```bash
python scripts/generate.py \
    --model-name evo-1_stripedhyena_pretrained_8k \
    --prompt ACGT \
    --n-samples 10 \
    --n-tokens 100 \
    --temperature 1. \
    --top-k 4 \
    --device cuda:0
```

We also provide an example script for using the model to score the log-likelihoods of a set of sequences.
```bash
python scripts/score.py \
    --input-fasta examples/example_seqs.fasta \
    --output-tsv scores.tsv \
    --model-name evo-1_stripedhyena_pretrained_8k \
    --device cuda:0
```
