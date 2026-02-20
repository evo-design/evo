import transformers
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM, BitsAndBytesConfig, AutoConfig, DataCollatorForLanguageModeling
import os

os.environ["CUDA_VISIBLE_DEVICES"] = "0"
os.environ["WANDB_DISABLED"] = "true"

model_name = 'togethercomputer/evo-1-8k-base'

model_config = AutoConfig.from_pretrained(model_name, trust_remote_code=True)
model_config.use_cache = False

model = AutoModelForCausalLM.from_pretrained(
    model_name,
    config=model_config,
    trust_remote_code=True,
    device_map={"":0},
    torch_dtype=torch.float16
)

tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
tokenizer.pad_token = "X"

for p in model.parameters():
    p.requires_grad = False

for p in model.backbone.blocks[-1].parameters():
    p.requires_grad = True

from datasets import load_dataset

dataset = load_dataset("gonzalobenegas/human-genome-cds")
data_collator = DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm=False)

def preprocess_function(sample):
    return tokenizer(sample['seq'], padding="longest", truncation=True, max_length=3000)

tokenized_ds = dataset.map(
    preprocess_function,
    batched=True,
    num_proc=12,
)

from transformers import AutoConfig, AutoModelForCausalLM, TrainingArguments, Trainer

training_args = TrainingArguments(
    output_dir="./evo_results",
    evaluation_strategy="epoch",
    learning_rate=2e-5,
    weight_decay=0.01,
    gradient_accumulation_steps=2,
    per_device_train_batch_size=4,
    warmup_steps=10,
    max_steps=100, # only a demo
    logging_steps=10,
    eval_steps=100,
    logging_strategy="steps",
    bf16=True
    # fp16=True, # This didn't work.
)


trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_ds["train"],
    eval_dataset=tokenized_ds["test"],
    data_collator=data_collator,
    
)

trainer.train()