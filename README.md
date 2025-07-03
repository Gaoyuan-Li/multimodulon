# MultiModulon

A Python package for analyzing multi-species/multi-strain/multi-modality profiles.

## Documentation

https://multimodulon.readthedocs.io/en/latest/

## Requirements

- Python >= 3.10
- NCBI BLAST+ tools (required for BBH generation)
- **PyTorch 2.6.0 with CUDA support (required for multi-view ICA)**
- **geotorch 0.3.0**

## Installation

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

# Create gene tables from GFF files
multiModulon.create_gene_table()

# Generate aligned expression matrices for multi-view ICA
multiModulon.generate_X('Output_Gene_Info')

# Print matrix dimensions and get recommendation for max components
# Shows aligned X matrices for all species with non-zero gene groups
# Provides maximum dimension recommendation for ICA

# STEP 1: Optimize the number of core components
optimal_num_core_components = multiModulon.optimize_number_of_core_components(
    metric='effect_size',          # Use Cohen's d effect size metric
    effect_size_threshold=5,       # Components must have Cohen's d > 5
    step=5,                        # Test k = 5, 10, 15, 20, ...
    save_path='optimization_plots', # Save plots to directory
    fig_size=(5, 3),              # Figure size
    font_path='/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'  # Optional custom font
)

print(f"Optimal number of core components: {optimal_num_core_components}")

# STEP 2: Optimize the number of unique components for each species
optimal_unique, optimal_total = multiModulon.optimize_number_of_unique_components(
    optimal_num_core_components=optimal_num_core_components,
    step=5,
    save_path='optimization_plots',  # Save all species plots
    fig_size=(5, 3),
    font_path='/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
)

print("Optimal unique components per species:", optimal_unique)
print("Optimal total components per species:", optimal_total)

# STEP 3: Run robust multi-view ICA with optimized parameters
M_matrices, A_matrices = multiModulon.run_robust_multiview_ica(
    a=optimal_total,           # Dictionary of total components per species
    c=optimal_num_core_components,  # Optimized core components
    num_runs=100,                   # Number of runs for robustness
    effect_size_threshold=5,        # Cohen's d threshold
    save_plots='ica_plots'          # Save clustering plots
)

# Alternative: Run standard multi-view ICA (single run)
multiModulon.run_multiview_ica(
    a=optimal_total,               # Can also use list or single value
    c=optimal_num_core_components  # Optimized core components
)

# Access the ICA results
for species_name in multiModulon._species_data.keys():
    M_matrix = multiModulon[species_name].M  # ICA mixing matrix
    A_matrix = multiModulon[species_name].A  # Activity matrix
    print(f"{species_name} - M: {M_matrix.shape}, A: {A_matrix.shape}")

# Calculate explained variance
explained_var = multiModulon.calculate_explained_variance()
for species, variance in explained_var.items():
    print(f"{species} explained variance: {variance:.4f}")

# Visualize iModulon weights
multiModulon.view_iModulon_weights(
    species='MG1655',
    component='Core_1',
    save_path='imodulon_plots',  # Save directory
    fig_size=(6, 4),
    font_path='/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
)

# Visualize core component weights across all species
multiModulon.view_core_iModulon_weights(
    component='Core_1',
    save_path='imodulon_plots',
    fig_size=(6, 4)
)

# Save the entire MultiModulon object
multiModulon.save_to_json_multimodulon('multimodulon_analysis.json')

# Load it back later
loaded_multiModulon = MultiModulon.load_json_multimodulon('multimodulon_analysis.json')
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

**Core Methods:**
- `__init__(input_folder_path)`: Initialize with Input_Data folder path
- `__getitem__(species_name)`: Get data for a specific strain
- `summary()`: Print data summary

**Data Preparation:**
- `generate_BBH(output_path, threads=1)`: Generate BBH files for all strain pairs
- `align_genes(input_bbh_dir, output_dir, reference_order, bbh_threshold)`: Align genes across strains
- `create_gene_table()`: Create gene tables from GFF files for all species
- `generate_X(gene_info_folder)`: Generate aligned expression matrices for multi-view ICA

**Optimization:**
- `optimize_number_of_core_components(metric='effect_size', save_path=None, fig_size=(5,3), font_path=None)`: 
  - Optimize number of core components using NRE or Cohen's d
  - Returns optimal number and saves plot as 'num_core_optimization.svg'
- `optimize_number_of_unique_components(optimal_num_core_components, save_path=None, fig_size=(5,3), font_path=None)`:
  - Optimize unique components per species
  - Saves plots as 'num_unique_{species}_optimization.svg'
  - Returns dictionaries of optimal unique and total components

**ICA Analysis:**
- `run_multiview_ica(a, c, mode='gpu')`: Run standard multi-view ICA (single run)
- `run_robust_multiview_ica(a, c, num_runs=100, save_plots=None)`: 
  - Run robust multi-view ICA with clustering
  - Returns M and A matrices
- `generate_A()`: Generate activity matrices from M matrices
- `calculate_explained_variance()`: Calculate explained variance for each species

**Visualization:**
- `view_iModulon_weights(species, component, save_path=None, fig_size=(6,4), font_path=None)`:
  - Visualize gene weights across genome for a component
  - Saves as '{species}_{component}_iModulon.svg'
- `view_core_iModulon_weights(component, save_path=None, fig_size=(6,4), font_path=None)`:
  - Visualize core component weights across all species
  - Saves individual plots for each species

**Data Management:**
- `get_orthologs(species1, species2)`: Get ortholog pairs (after BBH)
- `save_bbh(output_path)`: Save BBH results
- `load_bbh(input_path)`: Load BBH results
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

## Multi-view ICA Analysis

MultiModulon now supports GPU-accelerated multi-view Independent Component Analysis (ICA) for discovering shared and specific patterns across multiple species/strains.

### Key Features

- **Automatic Component Optimization**: Uses NRE (Noise Reduction Error) to automatically determine the optimal number of core components
- **GPU Acceleration**: Leverages PyTorch with CUDA support for fast computation
- **Flexible Input**: Works with any number of species/strains (≥2)
- **Orthogonal Constraints**: Uses geotorch for maintaining orthogonality in ICA unmixing matrices

### Complete Workflow

1. **Data Preparation**:
   - Generate BBH: `generate_BBH()` with multi-threading support
   - Align genes: `align_genes()` to create unified gene database
   - Create gene tables: `create_gene_table()` from GFF files
   - Generate aligned matrices: `generate_X()` for multi-view ICA

2. **Optimization**:
   - Find optimal core components: `optimize_number_of_core_components()`
   - Find optimal unique components: `optimize_number_of_unique_components()`

3. **ICA Analysis**:
   - Run robust ICA: `run_robust_multiview_ica()` with multiple runs
   - Or standard ICA: `run_multiview_ica()` for single run
   - A matrices are automatically generated from M matrices

4. **Analysis & Visualization**:
   - Calculate explained variance: `calculate_explained_variance()`
   - Visualize component weights: `view_iModulon_weights()` and `view_core_iModulon_weights()`
   - Save/load results: `save_to_json_multimodulon()` and `load_json_multimodulon()`

### Advanced Parameters

**Component Optimization:**
- `metric`: 'nre' or 'effect_size' (default: 'effect_size')
- `effect_size_threshold`: Cohen's d threshold (default: 5)
- `step`: Step size for testing k values (default: 5)
- `num_runs`: Number of cross-validation runs (default: 1)
- `train_frac`: Fraction of data for training (default: 0.75)
- `save_path`: Directory to save optimization plots
- `fig_size`: Figure size tuple (default: (5, 3))
- `font_path`: Path to custom font file

**Multi-view ICA:**
- `a`: Components per species (int, list, dict, or individual a1, a2, ...)
- `c`: Number of core components (use optimized value)
- `mode`: 'gpu' or 'cpu' computation mode (default: 'gpu')
- `effect_size_threshold`: For filtering components (default: 5)

**Robust Multi-view ICA:**
- `num_runs`: Number of ICA runs (default: 100)
- `effect_size_threshold_core`: Threshold for core components
- `effect_size_threshold_unique`: Threshold for unique components
- `num_top_gene`: Top genes for effect size calculation (default: 20)
- `save_plots`: Directory to save clustering plots

**Visualization:**
- `save_path`: Directory or file path for saving plots
- `fig_size`: Figure size as (width, height) in inches
- `font_path`: Path to TrueType font file for custom fonts

### Example Output

After running the complete workflow, you'll have:
- Aligned expression matrices (`species.X`) with consistent gene indexing
- ICA mixing matrices (`species.M`) revealing modular patterns
- Activity matrices (`species.A`) showing component activities
- Optimization plots:
  - Core components: `num_core_optimization.svg`
  - Unique components: `num_unique_{species}_optimization.svg`
- iModulon visualizations showing gene weights across genomes
- Explained variance values for each species
- Complete analysis saved in JSON format for reproducibility

## License

MIT License
