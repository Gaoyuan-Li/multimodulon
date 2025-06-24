# MultiModulon

A Python package for analyzing multi-species RNA-seq expression profiles with integrated ortholog detection.

## Features

- Load and manage expression data from multiple species
- Parse GFF annotations to extract gene information
- Perform Bidirectional Best Hits (BBH) analysis for ortholog detection
- Unified interface for accessing expression matrices and metadata
- Built-in data validation and sanity checks

## Installation

```bash
pip install -e .
```

For development:
```bash
pip install -e ".[dev]"
```

## Quick Start

```python
from multimodulon import MultiModulon

# Initialize with a data folder
multi_modulon = MultiModulon('/path/to/data/folder')

# Access species data
species_data = multi_modulon['Pseudomonas_fluorescens']

# Get expression matrices
log_tpm = species_data.log_tpm  # Raw log TPM values
X = species_data.X  # Normalized log TPM values

# Get sample metadata
sample_sheet = species_data.sample_sheet

# Get gene annotations
gene_table = species_data.gene_table

# Get orthologs between species
orthologs = multi_modulon.get_orthologs('Pseudomonas_fluorescens', 'Pseudomonas_putida')

# Save BBH results
multi_modulon.save_bbh('output/bbh_results')
```

## Data Structure

The package expects the following directory structure:

```
data_folder/
├── Species_name_1/
│   ├── expression_matrices/
│   │   ├── log_tpm.tsv
│   │   ├── log_tpm_norm.csv
│   │   ├── counts.tsv
│   │   └── tpm.tsv
│   ├── ref_genome/
│   │   ├── genomic.gff
│   │   └── *.fna
│   └── samplesheet/
│       └── samplesheet.csv
└── Species_name_2/
    └── ...
```

## BBH Analysis

The package uses BLAST+ via Docker for BBH analysis. Ensure Docker is installed and running:

```bash
docker pull quay.io/biocontainers/blast:2.16.0--h66d330f_5
```

## API Reference

### MultiModulon

Main class for multi-species analysis.

**Methods:**
- `__init__(data_folder, run_bbh=True)`: Initialize with data folder
- `__getitem__(species_name)`: Get data for a specific species
- `get_orthologs(species1, species2)`: Get ortholog pairs
- `save_bbh(output_path)`: Save BBH results
- `load_bbh(input_path)`: Load BBH results
- `summary()`: Print data summary

### SpeciesData

Container for single species data.

**Properties:**
- `log_tpm`: Log TPM expression matrix
- `X`: Normalized log TPM expression matrix
- `sample_sheet`: Sample metadata
- `gene_table`: Gene annotations

## License

MIT License