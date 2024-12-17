"""
Usage: 
python -m scripts.00_farm_activations --data_dir outputs --model evo-1-8k-base --batch_size 2 --max_seq_length 8192
"""
import torch
import numpy as np
import os
import json
import argparse
from datasets import load_dataset
from tqdm import tqdm
from pathlib import Path

from evo import Evo
from evo.scoring import prepare_batch

DATASET_ID = "LongSafari/open-genome"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, default='outputs',
                       help='Directory to save activation outputs')
    parser.add_argument('--model', type=str, default='evo-1-8k-base',
                       help='Model name/path to use')
    parser.add_argument('--batch_size', type=int, default=32,
                       help='Batch size for inference')
    parser.add_argument('--dataset_split', type=str, default='sample',
                       choices=['sample', 'stage1', 'stage2'],
                       help='Which split of open-genome to use')
    parser.add_argument('--max_seq_length', type=int, default=None,
                       help='Maximum sequence length to use')
    return parser.parse_args()

def register_hooks(model, data_dir, base_model):
    """Register forward hooks for each attention block."""
    hooks = []
    activations = []
    
    for name, module in model.named_modules():
        if 'AttentionBlock' in str(type(module)):
            # Create output directory for this module            
            def hook_fn(module, input, output, name=name):
                # Store the module name and save directory in the activation
                if isinstance(output, tuple):
                    activation = output[0]  # Usually the first element is the tensor we want
                else:
                    activation = output

                activations.append({
                    'name': name,
                    'output': activation.to(torch.float32).detach().cpu()
                })

            hooks.append(module.register_forward_hook(hook_fn))
    
    return hooks, activations

def save_batch_activations(activations, seq_lengths, data_dir, dataset_name, model_name, batch_idx):
    """Save activations for each module to separate files."""
    for act in activations:
        output_dir = Path(data_dir) / dataset_name / model_name / act['name'] 
        output_dir.mkdir(parents=True, exist_ok=True)

        output_path = output_dir / f"activations_{batch_idx}.pt"
        torch.save((act['output'], seq_lengths), output_path)

def save_metadata(args, data_dir, dataset_name, base_model):
    """Save metadata about the inference run."""
    metadata = {
        'batch_size': args.batch_size,
        'data_dir': data_dir,
        'base_model': base_model,
        'dataset_split': args.dataset_split
    }
    
    metadata_dir = Path(data_dir) / dataset_name / base_model
    metadata_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = metadata_dir / 'metadata.json'
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)

def run_inference(args):
    # Setup model
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    evo_model = Evo(args.model)
    model, tokenizer = evo_model.model, evo_model.tokenizer
    model.to(device)
    model.eval()

    # Setup data
    dataset_subsplit = "validation"
    dataset_name = f"{DATASET_ID}-{args.dataset_split}-{dataset_subsplit}".replace("/", "_")
    dataset = load_dataset(DATASET_ID, args.dataset_split)[dataset_subsplit]

    # Register hooks
    hooks, activations = register_hooks(model, args.data_dir, args.model)
    
    # Save metadata
    save_metadata(args, args.data_dir, dataset_name, args.model)
    
    # Run inference in batches
    try:
        for sample_idx in tqdm(range(0, len(dataset), args.batch_size)):
            batch = dataset[sample_idx:sample_idx + args.batch_size]
            batch_idx = sample_idx // args.batch_size
            
            # Prepare batch
            input_ids, seq_lengths = prepare_batch(
                batch['text'],
                tokenizer,
                max_seq_length=args.max_seq_length,
                prepend_bos=False,
                device=device
            )
            
            # Run inference and collect activations
            with torch.no_grad():
                _, _ = model(input_ids)  # Forward pass triggers hooks
                
                # Save activations for this batch
                save_batch_activations(activations, seq_lengths, args.data_dir, dataset_name, args.model, batch_idx)
                
                # Clear activations to free memory
                activations.clear()
                
    finally:
        # Clean up hooks
        for hook in hooks:
            hook.remove()

def main():
    args = parse_args()
    run_inference(args)

if __name__ == '__main__':
    main()    
