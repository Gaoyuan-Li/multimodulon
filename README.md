# MultiModulon

[![python](https://img.shields.io/badge/python-3.10_%7C_3.11_%7C_3.12-blue)](https://www.python.org)

A Python package for analyzing multi-species/multi-strain/multi-modality profiles.

[![PyTorch](https://img.shields.io/badge/PyTorch-EE4C2C?style=for-the-badge&logo=pytorch&logoColor=white)](https://pytorch.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-%23ffffff.svg?style=for-the-badge&logo=Matplotlib&logoColor=black)](https://matplotlib.org/stable/users/index.html)
[![Pandas](https://img.shields.io/badge/pandas-%23150458.svg?style=for-the-badge&logo=pandas&logoColor=white)](https://pandas.pydata.org/)
[![scikit-learn](https://img.shields.io/badge/scikit--learn-%23F7931E.svg?style=for-the-badge&logo=scikit-learn&logoColor=white)](https://scikit-learn.org/stable/)
[![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)](https://www.docker.com/)

## Documentation

https://multimodulon.readthedocs.io/en/latest/

## Requirements

- Python >= 3.10
- NCBI BLAST+ tools
- PyTorch (This package has been tested and verified to work with PyTorch 2.6.0 and 2.9.0)

## Installation

### Create multimodulon env

```bash
# Conda
conda create -n multimodulon python=3.11
conda activate multimodulon
```

### Install BLAST+ (Required for BBH generation)

```bash
# Conda
conda install -c bioconda blast
```

### Install PyTorch (Required for multi-view ICA)

```bash
# Install PyTorch with CUDA support
pip install torch==2.6.0 torchvision==0.21.0 torchaudio==2.6.0 --index-url https://download.pytorch.org/whl/cu124

# Install geotorch for orthogonal constraints
pip install geotorch==0.3.0
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

## Usage

Usage tutorial can be found at https://multimodulon.readthedocs.io/en/latest/examples/basic_workflow.html

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


## API Reference

### MultiModulon

Main class for multi-species/strain analysis.

**Constructor:**
- `__init__(input_folder_path)`: Initialize with Input_Data folder path

**Properties:**
- `bbh`: Get BBH (Bidirectional Best Hits) results
- `species`: Get list of loaded species names
- `combined_gene_db`: Get the combined gene database if available

**Data Access:**
- `__getitem__(species_name)`: Get data for a specific species (e.g., `multiModulon['species_name']`)
- `summary()`: Print summary of loaded data including species, samples, and matrix shapes

**Gene Analysis:**
- `generate_BBH(output_path="Output_BBH", threads=1)`: Generate BBH files using existing protein.faa files
- `align_genes(input_bbh_dir="Output_BBH", output_dir="Output_Gene_Info", reference_order=None, bbh_threshold=None)`: Align genes across all species using Union-Find algorithm
- `create_gene_table()`: Create gene tables from GFF files for all species
- `add_eggnog_annotation(eggnog_output_path)`: Add eggNOG annotations to gene tables
- `gff2pandas(gff_file, feature="CDS", index=None)`: Convert GFF file to pandas DataFrame

**Data Preparation:**
- `generate_X(gene_info_folder)`: Generate X matrices for all strains with consistent row indices
- `generate_A()`: Generate A matrices for all species from M matrices (A = M.T @ X)

**Optimization:**
- `optimize_number_of_core_components(**kwargs)`: Optimize number of core components using the single-gene filter metric
- `optimize_number_of_unique_components(**kwargs)`: Optimize unique components for each species
- `optimize_M_thresholds(method="Otsu's method", quantile_threshold=90)`: Optimize thresholds for M matrices

**ICA Analysis:**
- `run_multiview_ica(**kwargs)`: Run standard multi-view ICA (single run)
- `run_robust_multiview_ica(a, c, num_runs=100, ...)`: Run robust multi-view ICA with clustering
- `calculate_explained_variance()`: Calculate explained variance for each species

**iModulon Management:**
- `rename_iModulon(current_name, new_name, species=None)`: Rename an iModulon across all species and related dataframes

**Visualization:**
- `view_iModulon_weights(*args, **kwargs)`: View weights of an iModulon
- `view_iModulon_activities(*args, **kwargs)`: View activities of an iModulon
- `view_iModulon_genes(*args, **kwargs)`: View genes in an iModulon
- `view_core_iModulon_weights(*args, **kwargs)`: View weights of a core iModulon across all species
- `compare_core_iModulon(*args, **kwargs)`: Compare core iModulon across species
- `compare_core_iModulon_activity(*args, **kwargs)`: Compare core iModulon activities across species
- `plot_iM_conservation_bubble_matrix(*args, **kwargs)`: Summarize iModulon conservation across species
- `show_iModulon_activity_change(*args, **kwargs)`: Visualize iModulon activity changes between conditions
- `show_gene_iModulon_correlation(*args, **kwargs)`: Show correlation between gene expression and iModulon activity

**Data Management:**
- `get_orthologs(species1, species2)`: Get ortholog pairs between two species
- `save_bbh(output_path)`: Save BBH results to file
- `load_bbh(input_path)`: Load BBH results from file
- `save_to_json_multimodulon(save_path)`: Save entire MultiModulon object to JSON
- `load_json_multimodulon(load_path)`: Load MultiModulon object from JSON (static method)

### SpeciesData

Container for single strain data.

**Properties:**
- `log_tpm`: Log TPM expression matrix (from log_tpm.csv)
- `log_tpm_norm`: Normalized log TPM expression matrix (from log_tpm_norm.csv)
- `X`: Aligned expression matrix (used for multi-view ICA)
- `M`: ICA mixing matrix (available after running multi-view ICA)
- `A`: Activity matrix (components × samples, available after ICA)
- `sample_sheet`: Sample metadata
- `gene_table`: Gene annotations (parsed from GFF)


## License

MIT License
