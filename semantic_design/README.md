# Semantic Design

Semantic design harnesses genomic-context prompting with Evo to generate DNA or protein sequences that are enriched for targeted biological functions. The workflows in this directory reproduce the analyses described in the semantic-design manuscript, covering gene/operon completion as well as toxin–antitoxin and anti-CRISPR generation.

When run at scale, these pipelines can create (1) diversified variants of known proteins, (2) unannotated natural sequences that perform the same function, and (3) highly divergent candidates that retain functionality.

## Contents

```
./semantic_design/
├── bin/                                            # Miscellaneous helper scripts 
├── environments/                                   # Environments for running all pipelines
│   ├── semantic_design.yml                         # Environment for sampling and filtering with Evo
│   ├── esmfold_multimer.yml                        # Environment for running t2ta_cofold pipeline
├── pipelines/                                      # End-to-end sampling + filtering workflows
│   ├── gene_completion.py                          # Gene completion evaluation
│   ├── operon_completion.py                        # Operon completion evaluation
│   ├── t2ta_sample.py                              # Type II TA sampling + base filtering
│   ├── t2ta_cofold.py                              # ESMFold multimer + pDockQ scoring of TA pairs
│   ├── t3ta_sample.py                              # Type III TA sampling pipeline
│   └── acr_sample.py                               # Anti-CRISPR sampling pipeline
├── sample_configs/                                 # Minimal YAML configs for each workflow
├── sample_prompts/                                 # Example prompt CSVs referenced by sample configs
├── sample_reference_sequences/                     # Reference FASTA/CM/HMM files used in configs
├── semantic_design.py                              # Shared sampling/filtering utilities
├── semantic_design.yml                             # Conda environment for these workflows
└── README.md                                       # This document
```

## Setup

### Environment
1. Create and activate the conda environment (flash-attn install requires a GPU-enabled machine):
   ```bash
   conda env create -f semantic_design.yml
   conda activate semantic_design
   ```
2. Install Evo if it is not already available:
   ```bash
   pip install evo-model
   ```

### Running Pipelines
All scripts accept a YAML `--config` that specifies inputs and parameters (see `sample_configs/`). Each config requests a single `output_dir`; the pipeline will create the directory if needed and populate it with CSV/FASTA/MSA outputs.

Typical invocation:
```bash
conda activate semantic_design
python pipelines/<script>.py --config sample_configs/<config>.yaml
```

## Pipelines Overview

### `pipelines/gene_completion.py`
Evaluates Evo’s ability to complete partial sequences of highly conserved genes:
- Outputs gene-level summaries, per-completion FASTA files, and aggregate CSVs.
- Example config: `sample_configs/gene_completion.yaml`.

### `pipelines/operon_completion.py`
Tests multi-gene contextual design by completing operons:
- Generates MAFFT-guided alignments, filtered protein FASTAs, and summary statistics.
- Example config: `sample_configs/operon_completion.yaml`.

### `pipelines/t2ta_sample.py`
Semantic design for Type II toxin–antitoxin systems:
- Produces candidate CSV/FASTA files and Pfam-filtered protein sets ready for folding.
- Example config: `sample_configs/t2ta_sample.yaml`.

### `pipelines/t2ta_cofold.py`
Runs ESMFold multimer + pDockQ scoring on paired TA proteins emitted by the sampler.
- Requires activation of esmfold_multimer environment in environments subfolder.
- Example config: `sample_configs/t2ta_cofold.yaml`.

### `pipelines/t3ta_sample.py`
Semantic design for type III toxin-antitoxin generation:
- Generation → Prodigal/segmasker filtering → ESMFold triage.
- Tandem repeat finder, ViennaRNA folding, cmscan/hmmscan filters, and candidate TA pairing summaries.
- Example config: `sample_configs/t3ta_sample.yaml`.

### `pipelines/acr_sample.py`
Semantic design for anti-CRISPR discovery:
- Outputs candidate Acr sequences that can be scored with PaCRISPR or similar Acr prediction tools.
- Example config: `sample_configs/acr_sample.yaml`.

### `semantic_design.py`
Shared building blocks (prompt reading, Evo sampling, Prodigal/segmasker wrappers, FASTA writers) used by every pipeline. Import these utilities when creating new workflows.

## Example Data

- **Prompts**: `sample_prompts/` contains the CSVs referenced by each example config.
- **Reference sequences**: `sample_reference_sequences/` holds reference FASTA/CM/HMM files for alignment- or HMMScan and cmscan-based steps.

## Related Documentation

- [Main Evo README](https://github.com/evo-design/evo/blob/main/README.md) for installing Evo or running baseline generation/fine-tuning workflows.

## Citation

Please cite the semantic design paper when using these workflows:

```
@article{merchant2025semantic,
    author = {Merchant, Aditi T and King, Samuel H and Nguyen, Eric and Hie, Brian L},
    title = {Semantic design of functional de novo genes from a genomic language model},
    year = {2025},
    doi = {10.1038/s41586-025-09749-7},
    URL = {https://www.nature.com/articles/s41586-025-09749-7},
    journal = {Nature}
}
```
