# MultiModulon

A Python package for analyzing multi-species/multi-strain/multi-modality profiles.

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

# Generate aligned expression matrices for multi-view ICA
multiModulon.generate_X()

# Print matrix dimensions and get recommendation for max components
# Shows aligned X matrices for all species with non-zero gene groups
# Provides maximum dimension recommendation for ICA

# STEP 1: Optimize the number of core components using Cohen's d effect size
# This automatically finds the optimal number of core components
optimal_num_core_components, effect_scores = multiModulon.optimize_number_of_core_components(
    metric='effect_size',          # Use Cohen's d effect size metric
    effect_size_threshold=5,    # Components must have Cohen's d > 5
    step=5,                        # Test k = 5, 10, 15, 20, ...
    save_plot='optimization.png'   # Save optimization plot
)

print(f"Optimal number of core components: {optimal_num_core_components}")

# STEP 2: Run multi-view ICA with optimized parameters
# Multiple ways to specify component numbers:

# Option 1: Traditional way (for exactly 6 species)
if len(multiModulon._species_data) == 6:
    multiModulon.run_multiview_ica(
        a1=50, a2=50, a3=50, a4=50, a5=50, a6=50, 
        c=optimal_num_core_components  # Use optimized number of core components
    )

# Option 2: List format (for any number of species)
multiModulon.run_multiview_ica(
    a=[50, 50, 50],  # Components per species (adjust list length to match your species count)
    c=optimal_num_core_components  # Use optimized number of core components
)

# Option 3: Same number of components for all species
multiModulon.run_multiview_ica(
    a=50,                         # Same number of components for all species
    c=optimal_num_core_components  # Use optimized number of core components
)

# Access the ICA results (M matrices) for each species
for species_name in multiModulon._species_data.keys():
    M_matrix = multiModulon[species_name].M  # ICA mixing matrix
    print(f"{species_name} M matrix: {M_matrix.shape}")

# Access aligned expression matrices (same row indexes across all strains)
for species_name in multiModulon._species_data.keys():
    aligned_X = multiModulon[species_name].X  # Aligned expression matrix
    print(f"{species_name} X matrix: {aligned_X.shape}")
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
- **`generate_X()`**: Generate aligned expression matrices for multi-view ICA
- **`optimize_number_of_core_components(step=5, num_runs=3, train_frac=0.75, save_plot=None)`**: Optimize number of core components using NRE
- **`run_multiview_ica(**kwargs)`**: Run multi-view ICA on aligned expression matrices using PyTorch
- `get_orthologs(species1, species2)`: Get ortholog pairs (after BBH)
- `save_bbh(output_path)`: Save BBH results
- `load_bbh(input_path)`: Load BBH results
- `summary()`: Print data summary

### SpeciesData

Container for single strain data.

**Properties:**
- `log_tpm`: Log TPM expression matrix (from log_tpm.csv)
- `log_tpm_norm`: Normalized log TPM expression matrix (from log_tpm_norm.csv)
- `X`: Aligned expression matrix (used for multi-view ICA)
- **`M`: ICA mixing matrix (available after running multi-view ICA)**
- `sample_sheet`: Sample metadata
- `gene_table`: Gene annotations (parsed from GFF)

## Multi-view ICA Analysis

MultiModulon now supports GPU-accelerated multi-view Independent Component Analysis (ICA) for discovering shared and specific patterns across multiple species/strains.

### Key Features

- **Automatic Component Optimization**: Uses NRE (Noise Reduction Error) to automatically determine the optimal number of core components
- **GPU Acceleration**: Leverages PyTorch with CUDA support for fast computation
- **Flexible Input**: Works with any number of species/strains (≥2)
- **Orthogonal Constraints**: Uses geotorch for maintaining orthogonality in ICA unmixing matrices

### Workflow

1. **Prepare Data**: Run `generate_X()` to create aligned expression matrices
2. **Optimize Parameters**: Use `optimize_number_of_core_components()` to find optimal core components
3. **Run ICA**: Execute `run_multiview_ica()` with optimized parameters
4. **Access Results**: Get ICA mixing matrices from each species' `.M` property

### Advanced Parameters

**Component Optimization:**
- `step`: Step size for testing k values (default: 5)
- `num_runs`: Number of cross-validation runs for robustness (default: 3)
- `train_frac`: Fraction of data for training (default: 0.75)
- `mode`: 'gpu' or 'cpu' computation mode (default: 'gpu')

**Multi-view ICA:**
- `a`: Number of components per species (int, list, or individual a1, a2, ...)
- `c`: Number of core components (use optimized value)
- `mode`: 'gpu' or 'cpu' computation mode (default: 'gpu')

### Example Output

After running the complete workflow, you'll have:
- Aligned expression matrices (`species.X`) with consistent gene indexing
- ICA mixing matrices (`species.M`) revealing modular patterns
- Optimization plot showing NRE vs number of core components
- Recommendations for maximum dimensions based on data structure

## License

MIT License