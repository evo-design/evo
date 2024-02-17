import torch

from .generation import generate
from .models import load_checkpoint
from .scoring import prepare_batch


ckpt_path = '/checkpoint/etnguyen/7b_striped_120k/global_step63500'

if torch.cuda.is_available():
    device = 'cuda'
else:
    device = 'cpu'


def test_batched_logits_savanna(pytestconfig):
    """
    Tests that logits do not change with batch size.

    TODO: This test currently fails when it should not!
    """
    test_sequence = 'ACTGGA'
    
    if torch.cuda.is_available():
        device = 'cuda'
    else:
        device = 'cpu'

    model, tokenizer = load_checkpoint(ckpt_path, 'savanna')
    model = model.to(device)
    model = model.eval()
    
    with torch.inference_mode():
        input_ids_1x, _ = prepare_batch([test_sequence], tokenizer, device=device)
        logits_1x = model((input_ids_1x, None, None))
        if pytestconfig.getoption('verbose') > 0:
            print(logits_1x)

        input_ids_4x, _ = prepare_batch([test_sequence] * 4, tokenizer, device=device)
        logits_4x = model((input_ids_4x, None, None))
        if pytestconfig.getoption('verbose') > 0:
            print(logits_4x)

    assert torch.allclose(logits_1x[0], logits_4x[0]), \
        'Logits differ between different batch sizes.'


def test_compare_logits_savanna_recurrent_inference(pytestconfig):
    """
    Compare the logits between the two repos for performing inference.

    TODO: Unify inference code into a single repository.
    """
    test_sequence = 'ACTGGA'

    if torch.cuda.is_available():
        device = 'cuda'
    else:
        device = 'cpu'

    model_savanna, tokenizer_savanna = load_checkpoint(ckpt_path, 'savanna')
    (
        model_recurrent_inference,
        tokenizer_recurrent_inference,
    ) = load_checkpoint(ckpt_path, 'recurrent_inference')

    input_ids_savanna = torch.tensor(
        tokenizer_savanna.tokenize(test_sequence),
        dtype=torch.long,
    ).to(device).unsqueeze(0)
    input_ids_recurrent_inference = torch.tensor(
        tokenizer_recurrent_inference.tokenize(test_sequence),
        dtype=torch.long,
    ).to(device).unsqueeze(0)

    model_savanna.to(device)
    model_recurrent_inference.to(device)

    with torch.inference_mode():
        logits_savanna = model_savanna((input_ids_savanna, None, None))
        logits_recurrent_inference, _ = model_recurrent_inference(
            input_ids_recurrent_inference
        )

    if pytestconfig.getoption('verbose') > 0:
        print(logits_savanna)
        print(logits_recurrent_inference)

    assert torch.allclose(logits_savanna, logits_recurrent_inference), \
        'Logits differ between savanna and recurrent-inference.'


def test_batched_sampling(pytestconfig):
    """
    Various tests related to batched sampling.
    """
    model, tokenizer = load_checkpoint(ckpt_path, 'recurrent_inference')
    model = model.to(device)
    model = model.eval()

    # Test ability to generate without a prompt.
    assert generate(
        [ '', ],
        model,
        tokenizer,
        n_tokens=32,
        source_repository='recurrent_inference',
        temperature=0.,
        top_k=1,
        top_p=1.,
        verbose=1,
        device='cuda:0',
    ) is not None, 'Failed to sample from an empty prompt.'

    # Test ability to generate from prompts of different lengths.
    assert generate(
        [ '', 'A', 'AA', 'AAAA', ],
        model,
        tokenizer,
        n_tokens=32,
        source_repository='recurrent_inference',
        temperature=0.,
        top_k=1,
        top_p=1.,
        verbose=1,
        device='cuda:0',
    ) is not None, 'Failed to sample from different-length prompts.'

    generated_seqs, _ = generate(
        [ 'GCTCAA', 'GCTCAA','GCTCAA','GCTCAA' ],
        model,
        tokenizer,
        n_tokens=32,
        source_repository='recurrent_inference',
        temperature=0.,
        top_k=1,
        top_p=1.,
        verbose=1,
        device='cuda:0',
    )
    assert all(seq == generated_seqs[0] for seq in generated_seqs), \
        'Greedy0sampled different sequences for the same prompt.'

    generated_seqs_1x, _ = generate(
        [ 'GCTCAA', ],
        model,
        tokenizer,
        n_tokens=32,
        source_repository='recurrent_inference',
        temperature=0.,
        top_k=1,
        top_p=1.,
        verbose=1,
        device='cuda:0',
    )
    assert generated_seqs_1x[0] == generated_seqs[0], \
        'Different batch sizes led to different samples.'


#def test_compare_sampling_savanna_recurrent_inference(pytestconfig):
#    """
#    Compare the sampling procedures between the two repos for
#    performing inference.
#
#    TODO: Unify inference code into a single repository.
#    """
#
#    prompt_seqs = [ '', 'ATCG' ]
#
#    seqs_savanna = generate(prompt_seqs, 'savanna')
#
#    seqs_ri = generate(prompt_seqs, 'recurrent_inference')
