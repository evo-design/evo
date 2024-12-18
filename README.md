# Evo: DNA foundation modeling from molecular to genome scale

![Evo](evo.jpg)

Evo is a biological foundation model capable of long-context modeling and design.
Evo uses the [StripedHyena architecture](https://github.com/togethercomputer/stripedhyena) to enable modeling of sequences at a single-nucleotide, byte-level resolution with near-linear scaling of compute and memory relative to context length.
Evo has 7 billion parameters and is trained on [OpenGenome](https://huggingface.co/datasets/LongSafari/open-genome), a prokaryotic whole-genome dataset containing ~300 billion tokens.

We describe Evo in the paper [“Sequence modeling and design from molecular to genome scale with Evo”](https://www.science.org/doi/10.1126/science.ado9336).

We describe Evo 1.5 in the paper [“Semantic mining of functional _de novo_ genes from a genomic language model”](https://www.biorxiv.org/content/10.1101/2024.12.17.628962). We used the Evo 1.5 model to generate [SynGenome](https://evodesign.org/syngenome/), the first AI-generated genomics database containing over 100 billion base pairs of synthetic DNA sequences.

We provide the following model checkpoints:
| Checkpoint Name                        | Description |
|----------------------------------------|-------------|
| `evo-1.5-8k-base`   | A model pretrained with 8,192 context obtained by extending the pretraining of `evo-1-8k-base` to process 50% more training data. |
| `evo-1-8k-base`     | A model pretrained with 8,192 context. We use this model as the base model for molecular-scale finetuning tasks. |
| `evo-1-131k-base`   | A model pretrained with 131,072 context using `evo-1-8k-base` as the base model. We use this model to reason about and generate sequences at the genome scale. |
| `evo-1-8k-crispr`   | A model finetuned using `evo-1-8k-base` as the base model to generate CRISPR-Cas systems. |
| `evo-1-8k-transposon`   | A model finetuned using `evo-1-8k-base` as the base model to generate IS200/IS605 transposons. |

## News

**December 17, 2024:** We have found and fixed a bug in the code for Evo model inference affecting package versions from Nov 15-Dec 16, 2024, which has been corrected in release versions 0.3 and above. If you installed the package during this timeframe, please upgrade to correct the issue.

## Contents

- [Setup](#setup)
  - [Requirements](#requirements)
  - [Installation](#installation)
- [Usage](#usage)
- [HuggingFace](#huggingface)
- [Together API](#together-api)
- [colab](https://colab.research.google.com/github/evo-design/evo/blob/main/scripts/hello_evo.ipynb)
- [Playground wrapper](https://evo.nitro.bio/)
- [Dataset](#dataset)
- [Citation](#citation)

## Setup

### Requirements

Evo is based on [StripedHyena](https://github.com/togethercomputer/stripedhyena/tree/main).

Evo uses [FlashAttention-2](https://github.com/Dao-AILab/flash-attention), which may not work on all GPU architectures.
Please consult the [FlashAttention GitHub repository](https://github.com/Dao-AILab/flash-attention#installation-and-features) for the current list of supported GPUs.

Make sure to install the correct [PyTorch version](https://pytorch.org/) on your system.

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

We recommend that you install the PyTorch library first, before installing all other dependencies (due to dependency issues of the `flash-attn` library; see, e.g., this [issue](https://github.com/Dao-AILab/flash-attention/issues/246)).

One of our [example scripts](scripts/), demonstrating how to go from generating sequences with Evo to folding proteins ([scripts/generation_to_folding.py](scripts/generation_to_folding.py)), further requires the installation of `prodigal`. We have created an [environment.yml](environment.yml) file for this:

```bash
conda env create -f environment.yml
conda activate evo-design
```

## Usage

Below is an example of how to download Evo and use it locally through the Python API.
```python
from evo import Evo
import torch

device = 'cuda:0'

evo_model = Evo('evo-1-131k-base')
model, tokenizer = evo_model.model, evo_model.tokenizer
model.to(device)
model.eval()

sequence = 'ACGT'
input_ids = torch.tensor(
    tokenizer.tokenize(sequence),
    dtype=torch.int,
).to(device).unsqueeze(0)

with torch.no_grad():
    logits, _ = model(input_ids) # (batch, length, vocab)

print('Logits: ', logits)
print('Shape (batch, length, vocab): ', logits.shape)
```
An example of batched inference can be found in [`scripts/example_inference.py`](scripts/example_inference.py).

We provide an [example script](scripts/generate.py) for how to prompt the model and sample a set of sequences given the prompt.
```bash
python -m scripts.generate \
    --model-name 'evo-1-131k-base' \
    --prompt ACGT \
    --n-samples 10 \
    --n-tokens 100 \
    --temperature 1. \
    --top-k 4 \
    --device cuda:0
```

We also provide an [example script](scripts/score.py) for using the model to score the log-likelihoods of a set of sequences.
```bash
python -m scripts.score \
    --input-fasta examples/example_seqs.fasta \
    --output-tsv scores.tsv \
    --model-name 'evo-1-131k-base' \
    --device cuda:0
```

## HuggingFace

Evo is integrated with [HuggingFace](https://huggingface.co/togethercomputer/evo-1-131k-base).
```python
from transformers import AutoConfig, AutoModelForCausalLM

model_name = 'togethercomputer/evo-1-8k-base'

model_config = AutoConfig.from_pretrained(model_name, trust_remote_code=True, revision="1.1_fix")
model_config.use_cache = True

model = AutoModelForCausalLM.from_pretrained(
    model_name,
    config=model_config,
    trust_remote_code=True,
    revision="1.1_fix"
)
```


## Together API

Evo is available through Together AI with a [web UI](https://api.together.xyz/playground/language/togethercomputer/evo-1-131k-base), where you can generate DNA sequences with a chat-like interface.

For more detailed or batch workflows, you can call the Together API with a simple example below.


```python
import openai
import os

# Fill in your API information here.
client = openai.OpenAI(
  api_key=TOGETHER_API_KEY,
  base_url='https://api.together.xyz',
)

chat_completion = client.chat.completions.create(
  messages=[
    {
      "role": "system",
      "content": ""
    },
    {
      "role": "user",
      "content": "ACGT", # Prompt the model with a sequence.
    }
  ],
  model="togethercomputer/evo-1-131k-base",
  max_tokens=128, # Sample some number of new tokens.
  logprobs=True
)
print(
    chat_completion.choices[0].logprobs.token_logprobs,
    chat_completion.choices[0].message.content
)
```

## Dataset

The OpenGenome dataset for pretraining Evo is available at [Hugging Face datasets](https://huggingface.co/datasets/LongSafari/open-genome).

## Citation

Please cite the following publication when referencing Evo.

```
@article{nguyen2024sequence,
   author = {Eric Nguyen and Michael Poli and Matthew G. Durrant and Brian Kang and Dhruva Katrekar and David B. Li and Liam J. Bartie and Armin W. Thomas and Samuel H. King and Garyk Brixi and Jeremy Sullivan and Madelena Y. Ng and Ashley Lewis and Aaron Lou and Stefano Ermon and Stephen A. Baccus and Tina Hernandez-Boussard and Christopher Ré and Patrick D. Hsu and Brian L. Hie },
   title = {Sequence modeling and design from molecular to genome scale with Evo},
   journal = {Science},
   volume = {386},
   number = {6723},
   pages = {eado9336},
   year = {2024},
   doi = {10.1126/science.ado9336},
   URL = {https://www.science.org/doi/abs/10.1126/science.ado9336},
}
```

Please cite the following publication when referencing Evo 1.5.

```
@article {merchant2024semantic,
   author = {Merchant, Aditi T and King, Samuel H and Nguyen, Eric and Hie, Brian L},
   title = {Semantic mining of functional de novo genes from a genomic language model},
   year = {2024},
   doi = {10.1101/2024.12.17.628962},
   publisher = {Cold Spring Harbor Laboratory},
   URL = {https://www.biorxiv.org/content/early/2024/12/18/2024.12.17.628962},
   journal = {bioRxiv}
}
```
