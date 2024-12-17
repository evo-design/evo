"""
Usage: 
torchrun --nproc_per_node=4 scripts/00_farm_activations.py --data_dir outputs --model evo-1-8k-base --batch_size 2 --max_seq_length 8192
"""
import time
import torch
import numpy as np
import os
import json
import argparse
from datasets import load_dataset
from tqdm import tqdm
from pathlib import Path
import torch.distributed as dist
from torch.nn import DataParallel as DDP

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

def setup_distributed():
    """Initialize distributed environment."""
    if 'LOCAL_RANK' not in os.environ:
        os.environ['LOCAL_RANK'] = '0'
    
    # Initialize process group
    dist.init_process_group(backend='gloo')
    rank = dist.get_rank()
    world_size = dist.get_world_size()
    torch.cuda.set_device(int(os.environ["LOCAL_RANK"]))
    
    return rank, world_size

def save_batch_activations(activations, seq_lengths, experiment_metadata, data_dir, dataset_name, model_name, timestamp, rank, batch_idx):
    """Save activations for each module to separate files."""
    for act in activations:
        output_dir = Path(data_dir) / dataset_name / model_name / act['name'] / timestamp / f"rank_{rank}"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        is_metadata_saved = (output_dir / 'metadata.json').exists()
        if not is_metadata_saved:
            save_metadata(experiment_metadata, act, output_dir)

        output_path = output_dir / f"activations_{batch_idx}.pt"
        torch.save((act['output'], seq_lengths), output_path)

def save_metadata(experiment_metadata, activation, output_dir):
    """Save metadata about the inference run."""
    metadata_path = output_dir / 'metadata.json'
    experiment_metadata['module_name'] = activation['name']
    experiment_metadata['embed_dim'] = activation['output'].shape[-1]
    with open(metadata_path, 'w') as f:
        json.dump(experiment_metadata, f, indent=4)

def run_inference(args):
    # Setup distributed
    rank, world_size = setup_distributed()
    
    # Generate timestamp only on rank 0 and broadcast to all processes
    if rank == 0:
        timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    else:
        timestamp = None
    
    # Broadcast timestamp from rank 0 to all processes
    if world_size > 1:
        if rank == 0:
            timestamp_tensor = torch.tensor([ord(c) for c in timestamp], dtype=torch.long).cuda()
        else:
            timestamp_tensor = torch.zeros(19, dtype=torch.long).cuda()  # Length of timestamp string
        dist.broadcast(timestamp_tensor, src=0)
        if rank != 0:
            timestamp = ''.join([chr(i) for i in timestamp_tensor.cpu().numpy()])

    # Setup model
    device = f'cuda:{rank}'
    evo_model = Evo(args.model)
    model, tokenizer = evo_model.model, evo_model.tokenizer
    model.to(device)
    model = DDP(model, device_ids=[rank])
    model.eval()

    # Setup data
    dataset_subsplit = "validation"
    dataset_name = f"{DATASET_ID}-{args.dataset_split}-{dataset_subsplit}".replace("/", "_")
    dataset = load_dataset(DATASET_ID, args.dataset_split)[dataset_subsplit]
    
    # Register hooks
    hooks, activations = register_hooks(model.module, args.data_dir, args.model)
    experiment_metadata = {
        'data_dir': args.data_dir,
        'dataset_split': args.dataset_split,
        'dataset_subsplit': dataset_subsplit,
        'batch_size': args.batch_size,
        'effective_batch_size': args.batch_size * world_size,
        'max_seq_length': args.max_seq_length,
        'base_model': args.model,
        'timestamp': timestamp,
        'rank': rank,
        'world_size': world_size
    }
    
    # Run inference in batches
    try:
        for batch_idx in tqdm(range(len(dataset) // (world_size * args.batch_size) + 1)):
            start_idx = batch_idx * world_size * args.batch_size + rank * args.batch_size
            end_idx = min(start_idx + args.batch_size, len(dataset))
            
            if start_idx >= len(dataset):
                break

            batch = dataset[start_idx:end_idx]
            
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
                save_batch_activations(
                    activations, seq_lengths, experiment_metadata, 
                    args.data_dir, dataset_name, args.model, 
                    timestamp, rank, batch_idx
                )
                
                # Clear activations to free memory
                activations.clear()
                
    finally:
        # Clean up hooks
        for hook in hooks:
            hook.remove()
        
        # Clear any remaining activations
        activations.clear()
        
        # Clean up CUDA cache
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        
        # Move model to CPU before destruction
        model.to('cpu')
        
        # Clean up process group
        if dist.is_initialized():
            dist.destroy_process_group()

def main():
    args = parse_args()
    run_inference(args)

if __name__ == '__main__':
    main()    
