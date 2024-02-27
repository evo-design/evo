import argparse
import subprocess
import torch
import numpy as np
import biotite.structure.io as bsio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from transformers import (
    AutoConfig,
    AutoModelForCausalLM,
    AutoTokenizer,
    EsmForProteinFolding,
    set_seed
)

from stripedhyena.tokenizer import CharLevelTokenizer


def main():

    # Load command-line arguments.

    ap = argparse.ArgumentParser()

    # generation args:
    default_prompt = (
        "|d__Bacteria;"
        +"p__Pseudomonadota;"
        +"c__Gammaproteobacteria;"
        +"o__Enterobacterales;"
        +"f__Enterobacteriaceae;"
        +"g__Escherichia;"
        +"s__Escherichia|"
    )
    ap.add_argument('--prompt', type=str, default=default_prompt, help='Prompt for generation')
    ap.add_argument("--model-name", type=str, default="togethercomputer/evo-1-pretrained-131k", help='Hugging Face model name')
    ap.add_argument('--temperature', type=float, default=1.0, help='Temperature during sampling')
    ap.add_argument('--top-k', type=int, default=4, help='Top K during sampling')
    ap.add_argument('--top-p', type=float, default=1., help='Top P during sampling')
    ap.add_argument('--cached-generation', type=bool, default=True, help='Use KV caching during generation')
    ap.add_argument("--max-new-tokens", type=int, default=1024, help='Max new tokens during sampling')
    ap.add_argument("--repetition-penalty", type=float, default=1.0, help='Repetition penalty during sampling')
    ap.add_argument("--penalty-alpha", type=float, default=0.0, help='Penalty alpha during sampling')
    # output args:
    ap.add_argument("--sequence-fasta", type=str, default='sequence.fasta', help='Sequence fasta file')
    ap.add_argument("--proteins-fasta", type=str, default='proteins.fasta', help='Proteins fasta file')
    ap.add_argument("--structure-pdb", type=str, default='structure.pdb', help='Structure PDB file')
    # misc args:
    ap.add_argument('--device', type=str, default='cuda:0', help='Device for generation')
    ap.add_argument('--verbose', type=int, default=1, help='Verbosity level')
    ap.add_argument('--seed', type=int, default=12345, help='Random seed')

    args = ap.parse_args()

    # Set seed.

    torch.manual_seed(args.seed) # pytorch random seed
    np.random.seed(args.seed) # numpy random seed
    set_seed(args.seed) # huggingface random seed

    # Load model config.

    model_config = AutoConfig.from_pretrained(args.model_name, trust_remote_code=True)
    model_config.use_cache = True

    # Load model.

    print(f'Loading {args.model_name}...')
    model = AutoModelForCausalLM.from_pretrained(
        args.model_name,
        config=model_config,
        trust_remote_code=True,
    )
    model = model.to(args.device)
    model.backbone = model.backbone.to(torch.bfloat16)

    # Make character-level tokenizer.

    tokenizer = CharLevelTokenizer(vocab_size=512)

    # Encode prompt.

    print(f'Prompting {args.model_name} with: ', args.prompt)
    prompt_ids = torch.tensor(tokenizer.tokenize(args.prompt)).to(torch.long).to(args.device)

    # Generate.

    print('Generating...')
    gen_token_ids = model.generate(
        prompt_ids.unsqueeze(0), # add batch dimension
        max_new_tokens=args.max_new_tokens,
        temperature=args.temperature,
        repetition_penalty=args.repetition_penalty,
        top_k=args.top_k,
        top_p=args.top_p,
        penalty_alpha=args.penalty_alpha,
        do_sample=args.temperature is not None,
        eos_token_id=tokenizer.eos_id,
        pad_token_id=tokenizer.pad_id,
        use_cache=args.cached_generation,
    )
    
    # Decode.
    
    dna_seq = tokenizer.detokenize_batch(gen_token_ids)[0]
    print('Generated DNA sequence: ', dna_seq)
    
    # Saving generated sequence to fasta.

    dna_seq_record = SeqRecord(Seq(dna_seq), id="evo-dna", description="DNA sequence generated by Evo.")
    with open(args.sequence_fasta, "w") as output_handle:
        SeqIO.write(dna_seq_record, output_handle, "fasta")
    print('Saved DNA sequence to: ', args.sequence_fasta)

    # Predict genes from sequence.

    print('Predicting genes with prodigal...')
    cmd = f'prodigal -i {args.sequence_fasta} -a {args.proteins_fasta} -o genes.gbk -p meta'
    subprocess.run(cmd, shell=True)

    # Load ESMFold.

    print('Loading ESMFold...')
    esmfold = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
    esmfold = esmfold.to(args.device)
    esmfold.esm = esmfold.esm.half()

    # Load ESMFold tokenizer.

    esmfold_tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    
    # Fold proteins.

    print('Folding proteins with EMSFold...')
    for i, protein_record in enumerate(SeqIO.parse(args.proteins_fasta, "fasta")):
        protein_seq = str(protein_record.seq)[:-1] # remove stop codon
        print('Protein sequence: ', protein_seq)

        with torch.inference_mode():
            esmfold_in = esmfold_tokenizer([protein_seq], return_tensors="pt", add_special_tokens=False)
            esmfold_out = esmfold(**esmfold_in.to(args.device))
            esmfold_out_pdb = esmfold.output_to_pdb(esmfold_out)[0]

        with open(args.structure_pdb, "w") as f:
            f.write(esmfold_out_pdb)

        protein_struct = bsio.load_structure(args.structure_pdb, extra_fields=["b_factor"])
        print('Folded protein: ', protein_struct)


if __name__ == "__main__":
    main()