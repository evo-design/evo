# Semantic Design

This directory contains scripts for semantic design - an approach that harnesses genomic context-based prompting to generate DNA sequences enriched for targeted biological functions. Using Evo, these tools can enable function-guided design of proteins through a process analogous to genome mining via guilt-by-association. When done at a high-throughput scale, semantic design may produce sequences that contain 1) diversified variants of existing proteins, 2) unannotated natural proteins with the same functionality as the target system, and 3) highly divergent proteins with retained functionality.

The scripts in this directory provide examples of analyses conducted using semantic design:
- Gene completion evaluation for highly conserved genes
- Operon completion analysis for multi-gene sequences
- Generation and filtering of toxin-antitoxin sequences
- Generation and filtering of anti-CRISPR sequences

## Setup

### Environment Setup
1. Create the conda environment using the provided configuration. To install flash-attn, you may need to 
connect to a GPU. Installation may take up to 3 hours:
```bash
conda env create -f semantic_design.yml
conda activate semantic_design
```

2. Install Evo if it is not already installed. This should take no more than a few minutes:
```bash
pip install evo-model
```

## Scripts

### pipelines/gene_completion.py
Evaluates the model's ability to complete partial sequences of highly conserved genes:
- Produces CSV file with summary statistics for all gene completions
- Configuration specified via config file described in script
- See example configurations in `sample_configs/gene_completion.yaml`

### pipelines/operon_completion.py
Assesses in-context genomic design ability through operon completion:
- Generates CSV containing summary statistics for all operon completions
- Configuration specified via config file described in script
- See example configurations in `sample_configs/operon_completion.yaml`

### pipelines/toxin_antitoxin_sample.py
Implements semantic design for toxin-antitoxin system design:
- Produces FASTA file with toxin-antitoxin candidates which can be fed into ESMFold/AlphaFold2 to co-fold candidates
- Configuration specified via config file described in script
- See example configurations in `sample_configs/toxin_antitoxin_sample.yaml`

### pipelines/toxin_antitoxin_cofold.py
Runs ESMFold multimer and pDockQ scoring on paired toxin/antitoxin proteins:
- Consumes the paired CSV emitted by `toxin_antitoxin_sample.py`
- Executes folding inside the specified ESM_Fold repo and summarizes pDockQ metrics
- Requires running inside the environment that provides the ESMFold + PyTorch stack
- See example configurations in `sample_configs/toxin_antitoxin_cofold.yaml`

### pipelines/type_iii_ta_sample.py
Implements the Type III toxin-antitoxin sampling and analysis workflow:
- Follows the standard generation→filter→ESMFold path, then runs tandem-repeat and RNA structure analyses
- Outputs TRF tables, RNA fold predictions, and candidate TA pairing summaries
- Configuration specified via YAML (see `sample_configs/type_iii_ta_sample.yaml`)

### pipelines/acr_sample.py
Implements semantic design for toxin-antitoxin system design:
- Produces csv file with Acr candidates for input into PaCRISPR
- Configuration specified via config file described in script
- See example configurations in `sample_configs/acr_sample.yaml`

### semantic_design.py
Contains several functions that can be used for sampling and filtering generated sequences:
- Functions can be imported into separate scripts to generate new pipelines for sampling and filtering

## Usage

All scripts use YAML configuration files to specify inputs and parameters. Example configurations are provided in the `sample_configs/` directory. Each config now asks for a single `output_dir` instead of numerous file paths—the scripts will create the directory (if necessary) and populate it with consistently named CSV/FASTA outputs. Each script should take roughly 15 minutes to run using the example prompts.

To run any script, connect to a GPU and run:
```bash
conda activate semantic_design
python pipelines/script_name.py --config path/to/config.yaml
```

## Example Configurations and Prompts

Example configuration files are provided in the `sample_configs/` directory. 

Example prompts are provided in the `sample_prompts/` directory.

Reference sequences for the operon and gene completion scripts are provided in the `sample_reference_sequences` directory.

## Related Documentation

For more comprehensive documentation, please refer to:
- The [Main Repository README](https://github.com/evo-design/evo/blob/main/README.md) for setting up generation with Evo

## Citation

Please cite the following publication when referencing semantic design or Evo 1.5.

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
