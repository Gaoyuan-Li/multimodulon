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
- NCBI BLAST+ tools (required for BBH generation)

## Installation

### Install BLAST+ (Required for BBH generation)

```bash
# Conda
conda install -c bioconda blast
```

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
multiModulon = MultiModulon('/path/to/Input_Data')

# Access strain data
strain_data = multiModulon['MG1655']

# Get expression matrices
log_tpm = strain_data.log_tpm  # Raw log TPM values
log_tpm_norm = strain_data.log_tpm_norm  # Normalized log TPM values
X = strain_data.X  # Alias for log_tpm_norm

# Get sample metadata
sample_sheet = strain_data.sample_sheet

# Generate BBH files for all strain pairs
# Use threads parameter to speed up BLAST computation (default: 1)
multiModulon.generate_BBH('Output_BBH', threads=8)

# Align genes across all strains and create unified expression matrices
combined_gene_db = multiModulon.align_genes(
    input_bbh_dir='Output_BBH',
    output_dir='Output_Gene_Info',
    reference_order=['MG1655', 'BL21', 'C', 'Crooks', 'W', 'W3110'],  # optional
    bbh_threshold=90  # optional, minimum PID threshold
)

# Access aligned expression matrices (same row indexes across all strains)
for strain in multiModulon.species:
    aligned_X = multiModulon[strain].X  # Aligned expression matrix
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
│   │   ├── protein.faa        # Protein sequences (used for BBH)
│   │   └── *.fna              # Genome sequence
│   └── samplesheet/
│       └── samplesheet.csv     # Sample metadata
└── Strain_2/
    └── ...
```

## BBH Analysis

The package can generate BBH (Bidirectional Best Hits) files for ortholog detection using BLAST+. The `generate_BBH()` method will:

1. Use existing protein.faa files from each strain's ref_genome directory
2. Create mappings from protein IDs (NCBI_GP) to locus_tags using genomic.gff files
   - Standard strains: Uses `locus_tag` attribute from GFF files
   - W3110 strain: Extracts JW-numbers (e.g., JW4367) from Note field format "ECK0001:JW4367:b0001"
   - Warnings are printed when fallback extraction methods are used, including the strain name
3. Run all-vs-all BLAST comparisons between species
4. Save results as CSV files with locus_tag-to-locus_tag mappings

The BBH output files will contain the following columns:
- gene: Query gene locus_tag
- subject: Subject gene locus_tag
- PID: Percent identity
- alnLength: Alignment length
- mismatchCount: Number of mismatches
- gapOpenCount: Number of gap openings
- queryStart/queryEnd: Query alignment coordinates
- subjectStart/subjectEnd: Subject alignment coordinates
- eVal: E-value
- bitScore: Bit score
- gene_length: Query gene length
- COV: Coverage (alnLength/gene_length)
- BBH: Bidirectional best hit indicator ("<=>")

**Note:** BLAST+ must be installed before running BBH analysis (see Installation section above).

## API Reference

### MultiModulon

Main class for multi-species/strain analysis.

**Methods:**
- `__init__(input_folder_path)`: Initialize with Input_Data folder path
- `__getitem__(species_name)`: Get data for a specific strain
- `generate_BBH(output_path, threads=1)`: Generate BBH files for all strain pairs (threads: number of CPU threads for BLAST)
- `align_genes(input_bbh_dir, output_dir, reference_order, bbh_threshold)`: Align genes across strains and create unified matrices
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