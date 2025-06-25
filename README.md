# MultiModulon

A Python package for analyzing multi-species/multi-strain/multi-modality profiles.

## Features

- Load and manage expression data from multiple species/strains
- Parse GFF annotations to extract gene information using locus_tags
- Generate Bidirectional Best Hits (BBH) for ortholog detection
- Align genes across species using Union-Find algorithm
- Create unified expression matrices with consistent gene indexing
- Built-in data validation and sanity checks

## Requirements

- Python >= 3.10
- BLAST+ (for BBH generation)

## Installation

### Install from GitHub

```bash
# Install directly from GitHub (includes all dependencies)
pip install git+https://github.com/Gaoyuan-Li/multimodulon.git

# Or install a specific branch/tag
pip install git+https://github.com/Gaoyuan-Li/multimodulon.git@main
```

### Install for development

```bash
# Clone the repository
git clone https://github.com/Gaoyuan-Li/multimodulon.git
cd multimodulon

# Install in editable mode
pip install -e .
```

## Quick Start

```python
from multimodulon import MultiModulon

# Initialize with Input_Data folder path
multi_modulon = MultiModulon('/path/to/Input_Data')

# Access strain data
strain_data = multi_modulon['MG1655']

# Get expression matrices
log_tpm = strain_data.log_tpm  # Raw log TPM values
log_tpm_norm = strain_data.log_tpm_norm  # Normalized log TPM values
X = strain_data.X  # Alias for log_tpm_norm

# Get sample metadata
sample_sheet = strain_data.sample_sheet

# Generate BBH files for all strain pairs
multi_modulon.generate_BBH('Output_BBH')

# Align genes across all strains and create unified expression matrices
combined_gene_db = multi_modulon.align_genes('Output_Gene_Info')

# Access aligned expression matrices (same row indexes across all strains)
for strain in multi_modulon.species:
    aligned_X = multi_modulon[strain].X  # Aligned expression matrix
    print(f"{strain}: {aligned_X.shape}")
```

## Data Structure

The package expects the following directory structure:

```
Input_Data/
├── Strain_1/
│   ├── expression_matrices/
│   │   ├── log_tpm.csv        # Log-transformed TPM values
│   │   └── log_tpm_norm.csv   # Normalized log TPM values
│   ├── ref_genome/
│   │   ├── genomic.gff        # Gene annotations
│   │   └── *.fna              # Genome sequence
│   └── samplesheet/
│       └── samplesheet.csv     # Sample metadata
└── Strain_2/
    └── ...
```

## BBH Analysis

The package can generate BBH files using BLAST+. Make sure BLAST+ is installed:

```bash
# Install BLAST+ (Ubuntu/Debian)
sudo apt-get install ncbi-blast+

# Or use conda
conda install -c bioconda blast
```

## API Reference

### MultiModulon

Main class for multi-species/strain analysis.

**Methods:**
- `__init__(input_folder_path)`: Initialize with Input_Data folder path
- `__getitem__(species_name)`: Get data for a specific strain
- `generate_BBH(output_path)`: Generate BBH files for all strain pairs
- `align_genes(output_path)`: Align genes across strains and create unified matrices
- `get_orthologs(species1, species2)`: Get ortholog pairs (after BBH)
- `save_bbh(output_path)`: Save BBH results
- `load_bbh(input_path)`: Load BBH results
- `summary()`: Print data summary

### SpeciesData

Container for single strain data.

**Properties:**
- `log_tpm`: Log TPM expression matrix (from log_tpm.csv)
- `log_tpm_norm`: Normalized log TPM expression matrix (from log_tpm_norm.csv)
- `X`: Normalized log TPM expression matrix (alias for log_tpm_norm)
- `sample_sheet`: Sample metadata
- `gene_table`: Gene annotations (parsed from GFF)

## License

MIT License