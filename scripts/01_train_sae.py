"""
Usage:
python -m scripts.01_train_sae \
    --data_dir outputs/LongSafari_open-genome-sample-validation/evo-1-8k-base/blocks.8/2024-12-16_23-53-19 \
    --batch_size 32 \
    --learning_rate 1e-3 \
    --l1_coefficient 1e-3 \
    --epochs 1 \
    --output_dir outputs/sae_models
"""
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import os
from pathlib import Path
import argparse
from tqdm import tqdm
import glob
import time
import json

class SparseAutoencoder(nn.Module):
    def __init__(self, input_dim, hidden_dim, tied_weights=True):
        super().__init__()
        self.encoder = nn.Linear(input_dim, hidden_dim, bias=True)
        if tied_weights:
            self.decoder = nn.Linear(hidden_dim, input_dim, bias=True)
            self.decoder.weight = nn.Parameter(self.encoder.weight.t())
        else:
            self.decoder = nn.Linear(hidden_dim, input_dim, bias=True)
        
        # Initialize with scaled random weights
        nn.init.xavier_normal_(self.encoder.weight)
        if not tied_weights:
            nn.init.xavier_normal_(self.decoder.weight)

    def forward(self, x):
        # Flatten input if necessary
        original_shape = x.shape
        if len(original_shape) > 2:
            x = x.reshape(-1, original_shape[-1])
        
        # Encode and decode
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        
        # Reshape back if necessary
        if len(original_shape) > 2:
            decoded = decoded.reshape(original_shape)
        
        return decoded, encoded

class ActivationsDataset(Dataset):
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)
        
        # Load metadata from rank_0 (assuming all ranks have same metadata)
        metadata_file = next(self.data_dir.glob("rank_*/metadata.json"))
        with open(metadata_file, 'r') as f:
            self.metadata = json.load(f)
        
        # Collect files from all rank directories
        self.data_files = []
        for rank_dir in sorted(self.data_dir.glob("rank_*")):
            if rank_dir.is_dir():
                rank_files = sorted(list(rank_dir.glob("*.pt")))
                self.data_files.extend(rank_files)
        
        if not self.data_files:
            raise ValueError(f"No activation files found in {data_dir}")
        
        self.embed_dim = self.metadata['embed_dim']
        self.seq_length = self.metadata['max_seq_length']
        self.file_batch_size = self.metadata['batch_size']
        
        # Calculate total number of sequences
        self.total_sequences = len(self.data_files) * self.file_batch_size * self.seq_length
    
    def __len__(self):
        return self.total_sequences
    
    def __getitem__(self, idx):
        # Calculate which file and position within file
        file_idx = idx // (self.file_batch_size * self.seq_length)
        remaining_idx = idx % (self.file_batch_size * self.seq_length)
        batch_idx = remaining_idx // self.seq_length
        seq_idx = remaining_idx % self.seq_length

        # Load file if needed
        activations, _ = torch.load(self.data_files[file_idx])
        
        # Return single activation vector
        return activations[batch_idx, seq_idx]

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, required=True,
                       help='Directory containing activation files')
    parser.add_argument('--batch_size', type=int, default=32,
                       help='Number of files to load per batch')
    parser.add_argument('--learning_rate', type=float, default=1e-3,
                       help='Learning rate')
    parser.add_argument('--l1_coefficient', type=float, default=1e-3,
                       help='L1 sparsity coefficient')
    parser.add_argument('--epochs', type=int, default=10,
                       help='Number of epochs')
    parser.add_argument('--output_dir', type=str, default='sae_models',
                       help='Directory to save trained models')
    return parser.parse_args()

def train_epoch(model, dataloader, optimizer, l1_coefficient, device):
    model.train()
    total_loss = 0
    
    for batch in tqdm(dataloader, desc="Training"):
        batch = batch.to(device)
        optimizer.zero_grad()
        
        reconstructed, encoded = model(batch)
        
        # Reconstruction loss
        recon_loss = nn.functional.mse_loss(reconstructed, batch)
        
        # L1 sparsity loss
        l1_loss = l1_coefficient * encoded.abs().mean()
        
        # Total loss
        loss = recon_loss + l1_loss
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
    
    return total_loss / len(dataloader)

def main():
    args = parse_args()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Setup dataset and dataloader
    dataset = ActivationsDataset(args.data_dir)

    dataloader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=2,
        drop_last=False  # Allow partial batches
    )
    
    # Initialize model
    model = SparseAutoencoder(
        input_dim=dataset.embed_dim,
        hidden_dim=dataset.embed_dim
    ).to(device)
    
    # Setup optimizer
    optimizer = optim.Adam(model.parameters(), lr=args.learning_rate)

    input_metadata = dataset.metadata
    base_model = input_metadata['base_model']
    module_name = input_metadata['module_name']
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")

    final_output_dir = Path(args.output_dir) / base_model / module_name / timestamp
    final_output_dir.mkdir(parents=True, exist_ok=True)

    # Save metadata
    output_metadata = {
        'data_dir': args.data_dir,
        'batch_size': args.batch_size,
        'learning_rate': args.learning_rate,
        'l1_coefficient': args.l1_coefficient,
        'epochs': args.epochs,
        'model': base_model,
        'module_name': module_name,
        'timestamp': timestamp,
    }
    with open(final_output_dir / "metadata.json", 'w') as f:
        json.dump(output_metadata, f, indent=4)
    
    # Training loop
    for epoch in range(args.epochs):
        avg_loss = train_epoch(model, dataloader, optimizer, args.l1_coefficient, device)
        print(f"Epoch {epoch+1}/{args.epochs}, Loss: {avg_loss:.6f}")
    
    # Save model
    torch.save(model.state_dict(), final_output_dir / "sae.pt")

if __name__ == '__main__':
    main()
