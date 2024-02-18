import torch

from evo import Evo

if __name__ == '__main__':
    device = 'cuda:0'

    evo_model = Evo('evo-1_stripedhyena_pretrained_8k')
    model, tokenizer = evo_model.model, evo_model.tokenizer
    model.to(device)
    model.eval()

    # Example single-sequence inference.

    sequence = 'ACGT'
    input_ids = torch.tensor(
        tokenizer.tokenize(sequence),
        dtype=torch.int,
    ).to(device).unsqueeze(0)
    logits, _ = model(input_ids) # (batch, length, vocab)

    print('Logits: ', logits)
    print('Shape (batch, length, vocab): ', logits.shape)

    # Example batched inference.

    sequences = [
        'ACGT',
        'A',
        'AAAAACCCCCGGGGGTTTTT',
    ]

    from evo.scoring import prepare_batch

    input_ids, seq_lengths = prepare_batch(
        sequences,
        tokenizer,
        prepend_bos=False,
        device=device,
    )
    logits, _ = model(input_ids) # (batch, length, vocab)

    print('Batch logits: ', logits)
    print('Batch shape (batch, length, vocab): ', logits.shape)
    
