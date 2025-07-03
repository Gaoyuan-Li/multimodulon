"""Core MultiModulon class for multi-species expression analysis."""

from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import json
import pickle
import warnings
from typing import Dict, Optional, Tuple, List, Union

from .species_data import SpeciesData
from .gene_alignment import generate_BBH, align_genes
from .optimization import optimize_number_of_core_components, optimize_number_of_unique_components
from .gff_utils import gff2pandas, create_gene_table
from .utils.bbh import BBHAnalyzer
from .utils.fasta_utils import extract_protein_sequences
from .multiview_ica import run_multiview_ica
from .multiview_ica_optimization import calculate_cohens_d_effect_size
from sklearn.cluster import HDBSCAN
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch
import os

logger = logging.getLogger(__name__)

# COG category color mapping
COG_COLORS = {
    # INFORMATION STORAGE AND PROCESSING
    'Translation, ribosomal structure and biogenesis': 'black',
    'RNA processing and modification': 'lightblue',
    'Transcription': 'sandybrown',
    'Replication, recombination and repair': 'fuchsia',
    'Chromatin structure and dynamics': 'mediumpurple',
    
    # CELLULAR PROCESSES AND SIGNALING
    'Cell cycle control, cell division, chromosome partitioning': 'y',
    'Nuclear structure': 'goldenrod',
    'Defense mechanisms': 'lightgray',
    'Signal transduction mechanisms': 'lime',
    'Cell wall/membrane/envelope biogenesis': 'mediumvioletred',
    'Cell motility': 'orchid',
    'Cytoskeleton': 'darkviolet',
    'Extracellular structures': 'darkgoldenrod',
    'Intracellular trafficking, secretion, and vesicular transport': 'saddlebrown',
    'Posttranslational modification, protein turnover, chaperones': 'skyblue',
    'Post-translational modification, protein turnover, and chaperones': 'skyblue',  # Alternative spelling
    
    # METABOLISM
    'Energy production and conversion': 'lightgreen',
    'Carbohydrate transport and metabolism': 'pink',
    'Amino acid transport and metabolism': 'red',
    'Nucleotide transport and metabolism': 'c',
    'Coenzyme transport and metabolism': 'green',
    'Lipid transport and metabolism': 'turquoise',
    'Inorganic ion transport and metabolism': 'blue',
    'Secondary metabolites biosynthesis, transport and catabolism': 'dodgerblue',
    'Secondary metabolites biosynthesis, transport, and catabolism': 'dodgerblue',  # Alternative spelling
    
    # POORLY CHARACTERIZED
    'Function unknown': 'slategray',
    
    # No annotation
    'No COG annotation': 'lightskyblue'
}

# COG category letter codes
COG_LETTER_CODES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'S': 'Function unknown'
}


class MultiModulon:
    """
    Main class for multi-species expression analysis.
    
    This class provides a unified interface for accessing expression data,
    gene annotations, and ortholog relationships across multiple species.
    """
    
    def __init__(self, input_folder_path: str):
        """
        Initialize MultiModulon object.
        
        Parameters
        ----------
        input_folder_path : str
            Path to the Input_Data folder containing species/strain subfolders
        """
        print(f"\nInitializing MultiModulon...")
        
        self.input_folder_path = Path(input_folder_path)
        if not self.input_folder_path.exists():
            raise ValueError(f"Input folder not found: {input_folder_path}")
        
        # Initialize storage
        self._species_data: dict[str, SpeciesData] = {}
        self._bbh: dict[tuple[str, str], pd.DataFrame] | None = None
        self._combined_gene_db: pd.DataFrame | None = None
        
        # Load species data
        self._load_species_data()
    
    def _load_species_data(self):
        """Load data for all species in the input folder."""
        # Find all species directories
        species_dirs = [d for d in self.input_folder_path.iterdir() 
                       if d.is_dir() and not d.name.startswith('.')]
        
        if not species_dirs:
            raise ValueError(f"No species directories found in {self.input_folder_path}")
        
        logger.info(f"Found {len(species_dirs)} species directories")
        
        # Load each species
        print(f"\nLoading from {self.input_folder_path}:")
        print("=" * 60)
        
        for species_dir in species_dirs:
            species_name = species_dir.name
            logger.info(f"Loading data for {species_name}")
            
            try:
                species_data = SpeciesData(species_name, species_dir)
                # Load required data
                _ = species_data.log_tpm  # Trigger loading
                _ = species_data.log_tpm_norm  # Trigger loading
                _ = species_data.sample_sheet  # Trigger loading
                
                # Print species information
                print(f"\n{species_name}:")
                print(f"  - Number of genes: {species_data.log_tpm.shape[0]}")
                print(f"  - Number of samples: {species_data.log_tpm.shape[1]}")
                
                # Validate data consistency
                if species_data.validate_data():
                    self._species_data[species_name] = species_data
                    print(f"  - Data validation: PASSED")
                    logger.info(f"Successfully loaded {species_name}")
                else:
                    print(f"  - Data validation: FAILED")
                    logger.warning(f"Skipping {species_name} due to validation failure")
                    
            except Exception as e:
                print(f"  - Error: {str(e)}")
                logger.error(f"Failed to load {species_name}: {str(e)}")
        
        print("\n" + "=" * 60)
        print(f"Successfully loaded {len(self._species_data)} species/strains/modalities")
        print("=" * 60)
    
    def _run_bbh_analysis(self):
        """Run bidirectional best hits analysis."""
        logger.info("Running BBH analysis...")
        
        # Prepare protein FASTA files for each species
        species_fasta_data = {}
        
        for species_name, species_data in self._species_data.items():
            # Check if we need to extract protein sequences
            fasta_path = species_data.data_path / "ref_genome" / "proteins.faa"
            
            if not fasta_path.exists():
                logger.info(f"Extracting protein sequences for {species_name}")
                genome_fasta = species_data.data_path / "ref_genome" / list(species_data.data_path.glob("ref_genome/*.fna"))[0]
                gff_file = species_data.data_path / "ref_genome" / "genomic.gff"
                
                try:
                    extract_protein_sequences(genome_fasta, gff_file, fasta_path)
                except Exception as e:
                    logger.error(f"Failed to extract proteins for {species_name}: {str(e)}")
                    continue
            
            species_fasta_data[species_name] = {'fasta_path': fasta_path}
        
        # Run BBH analysis
        try:
            bbh_analyzer = BBHAnalyzer()
            self._bbh = bbh_analyzer.run_bbh_all_pairs(species_fasta_data)
            logger.info("BBH analysis completed")
        except Exception as e:
            logger.error(f"BBH analysis failed: {str(e)}")
            self._bbh = {}
    
    def __getitem__(self, species_name: str) -> SpeciesData:
        """
        Get data for a specific species.
        
        Parameters
        ----------
        species_name : str
            Name of the species
            
        Returns
        -------
        SpeciesData
            Data container for the species
        """
        if species_name not in self._species_data:
            raise KeyError(f"Species '{species_name}' not found. Available species: {list(self._species_data.keys())}")
        
        return self._species_data[species_name]
    
    @property
    def bbh(self) -> dict[tuple[str, str], pd.DataFrame] | None:
        """Get BBH results."""
        return self._bbh
    
    @property
    def species(self) -> list[str]:
        """Get list of loaded species."""
        return list(self._species_data.keys())
    
    @property
    def combined_gene_db(self) -> pd.DataFrame | None:
        """Get the combined gene database if available."""
        return getattr(self, '_combined_gene_db', None)
    
    def generate_A(self):
        """
        Generate A matrices for all species from M matrices and X matrices.
        
        A = M.T @ X
        
        Where:
        - M is the mixing matrix (genes x components)
        - X is the expression matrix (genes x samples)
        - A is the activity matrix (components x samples)
        
        The resulting A matrix has:
        - Row indices: component names from M columns
        - Column indices: sample names from X columns
        """
        print("\nGenerating A matrices...")
        
        for species_name in self._species_data.keys():
            species_data = self._species_data[species_name]
            
            # Check if M matrix exists
            if species_data._M is None:
                logger.warning(f"M matrix not found for {species_name}. Skipping A generation.")
                continue
            
            # Get M and X matrices
            M = species_data.M
            X = species_data.X
            
            # Calculate A = M.T @ X
            A = M.T @ X
            
            # A should have component names as rows and sample names as columns
            A.index = M.columns  # Component names
            A.columns = X.columns  # Sample names
            
            # Store A matrix in species data
            species_data._A = A
            
            print(f"✓ Generated A matrix for {species_name}: {A.shape}")
        
        print("\nA matrix generation completed!")
    
    def save_bbh(self, output_path: str):
        """
        Save BBH results to file.
        
        Parameters
        ----------
        output_path : str
            Path to save BBH results
        """
        if self._bbh is None:
            logger.warning("No BBH results to save")
            return
        
        output_path = Path(output_path)
        
        # Convert to serializable format
        bbh_data = {}
        for (sp1, sp2), df in self._bbh.items():
            key = f"{sp1}___{sp2}"
            bbh_data[key] = df.to_dict('records')
        
        # Save as JSON
        with open(output_path.with_suffix('.json'), 'w') as f:
            json.dump(bbh_data, f, indent=2)
        
        # Also save as pickle for easier loading
        with open(output_path.with_suffix('.pkl'), 'wb') as f:
            pickle.dump(self._bbh, f)
        
        logger.info(f"BBH results saved to {output_path}")
    
    def load_bbh(self, input_path: str):
        """
        Load BBH results from file.
        
        Parameters
        ----------
        input_path : str
            Path to load BBH results from
        """
        input_path = Path(input_path)
        
        if input_path.with_suffix('.pkl').exists():
            # Load from pickle
            with open(input_path.with_suffix('.pkl'), 'rb') as f:
                self._bbh = pickle.load(f)
        elif input_path.with_suffix('.json').exists():
            # Load from JSON
            with open(input_path.with_suffix('.json'), 'r') as f:
                bbh_data = json.load(f)
            
            self._bbh = {}
            for key, records in bbh_data.items():
                sp1, sp2 = key.split('___')
                self._bbh[(sp1, sp2)] = pd.DataFrame(records)
        else:
            raise FileNotFoundError(f"BBH file not found: {input_path}")
        
        logger.info(f"BBH results loaded from {input_path}")
    
    def get_orthologs(self, species1: str, species2: str) -> pd.DataFrame:
        """
        Get ortholog pairs between two species.
        
        Parameters
        ----------
        species1, species2 : str
            Species names
            
        Returns
        -------
        pd.DataFrame
            Ortholog pairs with columns: query_id, subject_id, evalue, bitscore
        """
        if self._bbh is None:
            raise ValueError("BBH analysis not run. Initialize with run_bbh=True or load BBH results.")
        
        key = (species1, species2)
        if key not in self._bbh:
            raise ValueError(f"No BBH results found for {species1} vs {species2}")
        
        return self._bbh[key]
    
    def generate_X(self, gene_info_folder: str):
        """
        Generate X matrices for all strains with consistent row indices based on combined_gene_db.
        
        This method loads the combined gene database from the specified folder and creates
        aligned expression matrices (X) for all strains. All X matrices will have the same
        row indices, using gene names from the leftmost non-null entry in each gene group.
        The original log_tpm_norm data is preserved unchanged; only the X matrices are created
        with aligned gene names.
        
        Parameters
        ----------
        gene_info_folder : str
            Path to the Gene_Info folder containing combined_gene_db.csv
            Example: "../imminer_2_industrial_strain/Output_Gene_Info"
        """
        gene_info_path = Path(gene_info_folder)
        combined_db_path = gene_info_path / "combined_gene_db.csv"
        
        if not combined_db_path.exists():
            raise FileNotFoundError(f"combined_gene_db.csv not found in {gene_info_folder}")
        
        # Load the combined gene database
        logger.info(f"Loading combined gene database from {combined_db_path}")
        combined_gene_db = pd.read_csv(combined_db_path)
        self._combined_gene_db = combined_gene_db
        
        # Get list of strains from the columns
        strains = list(combined_gene_db.columns)
        logger.info(f"Found {len(strains)} strains in combined gene database: {strains}")
        
        # Verify all strains exist in loaded species data
        available_strains = set(self._species_data.keys())
        missing_strains = set(strains) - available_strains
        if missing_strains:
            logger.warning(f"The following strains in combined_gene_db are not loaded in MultiModulon: {missing_strains}")
            # Filter to only use available strains
            strains = [s for s in strains if s in available_strains]
            logger.info(f"Proceeding with {len(strains)} available strains: {strains}")
        
        # Extract gene names from the leftmost non-null entry in each row
        logger.info("Extracting gene names from leftmost entries in gene groups...")
        gene_names = []
        for idx, row in combined_gene_db.iterrows():
            # Find the first non-null gene in the row (from left to right)
            gene_name = None
            for strain in combined_gene_db.columns:
                val = row[strain]
                # Check for non-null and non-"None" values
                if pd.notna(val) and val != "None" and val is not None:
                    gene_name = val
                    break
            
            if gene_name is not None:
                gene_names.append(gene_name)
            else:
                # If all entries are null, use a placeholder
                gene_names.append(f"gene_group_{idx}")
        
        # Create a mapping from row index to gene name
        gene_index = pd.Index(gene_names)
        logger.info(f"Created gene index with {len(gene_index)} genes")
        
        # Create aligned X matrices for each strain
        logger.info("Creating aligned X matrices for all strains...")
        for strain in strains:
            if strain not in self._species_data:
                logger.warning(f"Strain {strain} not found in loaded species data, skipping")
                continue
            
            species_data = self._species_data[strain]
            
            # Get original expression matrix (preserve the original log_tpm_norm)
            orig_expr = species_data.log_tpm_norm
            if orig_expr.empty:
                logger.warning(f"No expression data found for {strain}, skipping")
                continue
            
            cols = orig_expr.columns
            
            # Create new aligned matrix with gene_index as row index
            # Only include genes that exist in combined_gene_db (no new zero rows)
            aligned_X = pd.DataFrame(index=gene_index, columns=cols, dtype=float)
            
            # Create a mapping from strain genes to leftmost gene names
            strain_gene_to_leftmost = {}
            for gene_name, row in zip(gene_names, combined_gene_db.itertuples(index=False)):
                # Get the gene ID for this strain in this gene group
                strain_gene_id = getattr(row, strain, None)
                if pd.notna(strain_gene_id) and strain_gene_id != "None" and strain_gene_id is not None:
                    strain_gene_to_leftmost[strain_gene_id] = gene_name
            
            # Fill in expression values by mapping original genes to leftmost gene names
            for gene_id in orig_expr.index:
                if gene_id in strain_gene_to_leftmost:
                    # This gene has a mapping to a gene group - use the leftmost name
                    leftmost_name = strain_gene_to_leftmost[gene_id]
                    aligned_X.loc[leftmost_name] = orig_expr.loc[gene_id]
            
            # Fill remaining rows with zeros (genes not present in this strain)
            aligned_X = aligned_X.fillna(0.0)
            
            # Update only the X matrix (preserve original log_tpm_norm)
            self._species_data[strain]._X = aligned_X
            
            # Log statistics
            non_zero_genes = (aligned_X != 0).any(axis=1).sum()
            logger.info(f"Created aligned X matrix for {strain}: shape={aligned_X.shape}, "
                       f"non-zero genes={non_zero_genes}/{len(aligned_X)}")
        
        logger.info("Successfully generated aligned X matrices for all strains")
        
        # Calculate maximum dimension recommendation
        min_columns = min([self._species_data[strain].X.shape[1] for strain in strains if strain in self._species_data])
        # Use the n*10 or n*10+5 that is smaller than the lowest number of columns
        max_dim = ((min_columns // 10) * 10) if min_columns % 10 < 5 else ((min_columns // 10) * 10 + 5)
        if max_dim >= min_columns:
            max_dim = ((min_columns // 10) * 10) if min_columns >= 10 else min_columns - 5
        
        # Print summary
        print(f"\nGenerated aligned X matrices:")
        print("=" * 60)
        for strain in strains:
            if strain in self._species_data:
                X = self._species_data[strain].X
                non_zero = (X != 0).any(axis=1).sum()
                print(f"{strain}: {X.shape} ({non_zero} non-zero gene groups)")
        print("=" * 60)
        print(f"Maximum dimension recommendation: {max_dim}")
        print("=" * 60)
    
    def summary(self):
        """Print summary of loaded data."""
        print(f"MultiModulon Summary")
        print(f"====================")
        print(f"Input folder: {self.input_folder_path}")
        print(f"Number of species/strains: {len(self._species_data)}")
        print(f"\nSpecies/strains loaded:")
        
        for species_name, species_data in self._species_data.items():
            print(f"\n  {species_name}:")
            try:
                print(f"    - Samples: {species_data.log_tpm.shape[1]}")
                print(f"    - Log TPM matrix shape: {species_data.log_tpm.shape}")
                print(f"    - Log TPM normalized matrix shape: {species_data.log_tpm_norm.shape}")
                print(f"    - X matrix shape (alias for log_tpm_norm): {species_data.X.shape}")
            except Exception as e:
                print(f"    - Error accessing data: {str(e)}")
    
    def run_multiview_ica(self, **kwargs):
        """
        Run multi-view ICA on aligned expression matrices.
        
        Parameters
        ----------
        a : int or list
            Number of components per species (int for all, list for each)
        c : int
            Number of core components
        mode : str, optional
            'gpu' or 'cpu' (default: 'gpu')
        effect_size_threshold : float, optional
            Cohen's d threshold for component filtering
        
        Examples
        --------
        >>> multiModulon.run_multiview_ica(a=50, c=30)  # same a for all
        >>> multiModulon.run_multiview_ica(a=[50, 50, 50], c=30)  # different a per species
        """
        # Check if X matrices have been generated
        species_list = list(self._species_data.keys())
        n_species = len(species_list)
        if n_species < 2:
            raise ValueError(f"Multi-view ICA requires at least 2 species/strains, found {n_species}")
        
        print(f"Running multi-view ICA on {n_species} species/strains: {species_list}")
        
        # Check if all species have X matrices
        for species in species_list:
            if self._species_data[species]._X is None:
                raise ValueError(
                    f"X matrix not found for {species}. "
                    "Please run generate_X() first to create aligned expression matrices."
                )
        
        # Extract a_values and c from kwargs
        a_values = {}
        c = kwargs.get('c')
        if c is None:
            raise ValueError("Parameter 'c' (number of core components) is required")
        
        # Extract a values for each species - support flexible naming
        for i, species in enumerate(species_list, 1):
            a_key = f'a{i}'
            if a_key in kwargs:
                a_values[species] = kwargs[a_key]
            elif 'a' in kwargs and isinstance(kwargs['a'], (list, tuple)):
                # Support passing a list of a values
                if len(kwargs['a']) != n_species:
                    raise ValueError(f"Length of 'a' list ({len(kwargs['a'])}) must match number of species ({n_species})")
                a_values[species] = kwargs['a'][i-1]
            elif 'a' in kwargs and isinstance(kwargs['a'], int):
                # Support passing a single a value for all species
                a_values[species] = kwargs['a']
            else:
                raise ValueError(f"Parameter '{a_key}' is required for species {species}, or provide 'a' as list/int")
        
        # Get other parameters
        mode = kwargs.get('mode', 'gpu')
        effect_size_threshold = kwargs.get('effect_size_threshold', None)
        
        # Prepare X matrices dictionary
        species_X_matrices = {}
        for species in species_list:
            species_X_matrices[species] = self._species_data[species].X
        
        # Run multi-view ICA
        results = run_multiview_ica(
            species_X_matrices=species_X_matrices,
            a_values=a_values,
            c=c,
            mode=mode,
            effect_size_threshold=effect_size_threshold
        )
        
        # Save M matrices to each species
        print("\nSaving M matrices to species objects...")
        for species, M_matrix in results.items():
            self._species_data[species].M = M_matrix
            print(f"✓ Saved M matrix for {species}: {M_matrix.shape}")
        
        # Generate A matrices from M matrices
        print("\nGenerating A matrices from M matrices...")
        for species in species_list:
            M = self._species_data[species].M
            X = self._species_data[species].X
            
            # Calculate A = M.T @ X
            A = M.T @ X
            
            # A should have component names as rows and sample names as columns
            A.index = M.columns  # Component names
            A.columns = X.columns  # Sample names
            
            # Store A matrix in species data
            self._species_data[species]._A = A
            
            print(f"✓ Generated A matrix for {species}: {A.shape}")
        
        print("\nMulti-view ICA completed successfully!")
    
    def run_robust_multiview_ica(
        self, 
        a: Dict[str, int],
        c: int,
        num_runs: int = 100,
        mode: str = 'gpu',
        seed: int = 42,
        effect_size_threshold: Optional[float] = 5,
        effect_size_threshold_core: Optional[float] = None,
        effect_size_threshold_unique: Optional[float] = None,
        num_top_gene: int = 20
    ) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
        """
        Run robust multi-view ICA with clustering to identify consistent components.
        
        This method runs multi-view ICA multiple times and uses clustering to identify
        robust components that appear consistently across runs. Components are filtered
        by Cohen's d effect size before clustering.
        
        Parameters
        ----------
        a : dict
            Dictionary mapping species names to total number of components
            Example: {'strain1': 50, 'strain2': 55, 'strain3': 60}
        c : int
            Number of core (shared) components
        num_runs : int, default=100
            Number of ICA runs to perform
        mode : str, default='gpu'
            'gpu' or 'cpu' mode for computation
        seed : int, default=42
            Random seed for reproducibility
        effect_size_threshold : float, optional, default=5
            Cohen's d threshold for both core and unique components.
            Only used if specific thresholds are not provided.
        effect_size_threshold_core : float, optional
            Cohen's d threshold specifically for core components.
            If provided, overrides effect_size_threshold for core components.
        effect_size_threshold_unique : float, optional
            Cohen's d threshold specifically for unique components.
            If provided, overrides effect_size_threshold for unique components.
        num_top_gene : int, default=20
            Number of top genes to use when calculating Cohen's d effect size
            
        Returns
        -------
        M_matrices : dict
            Dictionary mapping species names to M matrices containing robust components.
            Each M matrix has columns for core components followed by unique components.
        A_matrices : dict
            Dictionary mapping species names to A matrices (component activities).
            Each A matrix has components as rows and samples as columns.
            
        Examples
        --------
        >>> # Run with optimized parameters
        >>> a_values = {'strain1': 50, 'strain2': 55, 'strain3': 60}
        >>> M_matrices, A_matrices = multiModulon.run_robust_multiview_ica(
        ...     a=a_values,
        ...     c=30,
        ...     num_runs=100,
        ...     effect_size_threshold=5
        ... )
        
        >>> # Use different thresholds for core and unique components
        >>> M_matrices, A_matrices = multiModulon.run_robust_multiview_ica(
        ...     a=a_values,
        ...     c=30,
        ...     effect_size_threshold_core=5,
        ...     effect_size_threshold_unique=3
        ... )
        
        >>> # Access the results
        >>> M_strain1 = M_matrices['strain1']  # Mixing matrix for strain1
        >>> A_strain1 = A_matrices['strain1']  # Activity matrix for strain1
        """
        # Set effect size thresholds
        if effect_size_threshold_core is None:
            effect_size_threshold_core = effect_size_threshold if effect_size_threshold is not None else 5
        if effect_size_threshold_unique is None:
            effect_size_threshold_unique = effect_size_threshold if effect_size_threshold is not None else 5
            
        species_list = list(self._species_data.keys())
        n_species = len(species_list)
        
        # Validate inputs
        if n_species < 2:
            raise ValueError(f"Multi-view ICA requires at least 2 species/strains, found {n_species}")
        
        if not all(species in a for species in species_list):
            missing = [s for s in species_list if s not in a]
            raise ValueError(f"Missing 'a' values for species: {missing}")
            
        print(f"\nRunning robust multi-view ICA with {num_runs} runs")
        print(f"Species: {species_list}")
        print(f"Total components (a): {a}")
        print(f"Core components (c): {c}")
        print(f"Effect size thresholds - Core: {effect_size_threshold_core}, Unique: {effect_size_threshold_unique}")
        
        # Check if all species have X matrices
        for species in species_list:
            if self._species_data[species]._X is None:
                raise ValueError(
                    f"X matrix not found for {species}. "
                    "Please run generate_X() first to create aligned expression matrices."
                )
        
        # Prepare X matrices
        species_X_matrices = {}
        for species in species_list:
            species_X_matrices[species] = self._species_data[species].X
        
        # Initialize storage for components
        core_components = {species: [] for species in species_list}
        unique_components = {species: [] for species in species_list}
        
        # Run ICA multiple times
        print(f"\nCollecting components from {num_runs} runs...")
        for run_idx in tqdm(range(num_runs), desc="ICA runs"):
            # Run multi-view ICA with current seed
            M_matrices = run_multiview_ica(
                species_X_matrices=species_X_matrices,
                a_values=a,
                c=c,
                mode=mode,
                seed=seed + run_idx  # Different seed for each run
            )
            
            # Process components for each species
            for species in species_list:
                M = M_matrices[species]
                
                # Process core components
                for comp_idx in range(c):
                    weight_vector = M.iloc[:, comp_idx].values
                    
                    # Calculate effect size
                    effect_size = calculate_cohens_d_effect_size(weight_vector, seed, num_top_gene)
                    
                    # Only keep components above threshold
                    if effect_size >= effect_size_threshold_core:
                        # Apply sign convention
                        component = self._enforce_sign_convention(weight_vector)
                        core_components[species].append(component)
                
                # Process unique components  
                for comp_idx in range(c, a[species]):
                    weight_vector = M.iloc[:, comp_idx].values
                    
                    # Calculate effect size
                    effect_size = calculate_cohens_d_effect_size(weight_vector, seed, num_top_gene)
                    
                    # Only keep components above threshold
                    if effect_size >= effect_size_threshold_unique:
                        # Apply sign convention
                        component = self._enforce_sign_convention(weight_vector)
                        unique_components[species].append(component)
        
        print(f"\nClustering components...")
        
        # Cluster core components across all species
        core_results = self._cluster_core_components(
            core_components, species_list, num_runs, n_species
        )
        
        # Cluster unique components per species
        unique_results = self._cluster_unique_components(
            unique_components, species_list, num_runs
        )
        
        # Create final M matrices with robust components
        print(f"\nCreating final M matrices with robust components...")
        final_M_matrices = {}
        
        for species in species_list:
            # Get gene names from X matrix
            gene_names = species_X_matrices[species].index
            
            # Collect all components for this species
            components = []
            component_names = []
            
            # Add core components
            for core_name, centroids in core_results['centroids'].items():
                if species in centroids:
                    components.append(centroids[species])
                    component_names.append(core_name)
            
            # Add unique components
            for unique_name, centroid in unique_results[species]:
                components.append(centroid)
                component_names.append(unique_name)
            
            if components:
                # Create M matrix
                M = pd.DataFrame(
                    np.column_stack(components),
                    index=gene_names,
                    columns=component_names
                )
                final_M_matrices[species] = M
            else:
                # Create empty M matrix if no components found
                final_M_matrices[species] = pd.DataFrame(
                    index=gene_names,
                    columns=[]
                )
        
        # Save M matrices to species objects
        print("\nSaving robust M matrices to species objects...")
        for species, M_matrix in final_M_matrices.items():
            self._species_data[species].M = M_matrix
            n_core = len([c for c in M_matrix.columns if c.startswith('Core_')])
            n_unique = len([c for c in M_matrix.columns if c.startswith('Unique_')])
            print(f"✓ {species}: {M_matrix.shape} ({n_core} core, {n_unique} unique components)")
        
        # Generate A matrices using the robust M matrices
        print("\nGenerating A matrices from robust M matrices...")
        final_A_matrices = {}
        
        for species in species_list:
            M = final_M_matrices[species]
            X = species_X_matrices[species]
            
            # Calculate A = M.T @ X
            A = M.T @ X
            
            # A should have component names as rows and sample names as columns
            A.index = M.columns  # Component names
            A.columns = X.columns  # Sample names
            
            # Store A matrix
            final_A_matrices[species] = A
            self._species_data[species]._A = A
            
            print(f"✓ Generated A matrix for {species}: {A.shape}")
        
        # Print summary
        print(f"\n{'='*60}")
        print("Robust multi-view ICA completed!")
        print(f"{'='*60}")
        print(f"Total core clusters identified: {len(core_results['clusters'])}")
        for species in species_list:
            n_unique = len(unique_results[species])
            print(f"{species}: {n_unique} unique components")
        print(f"{'='*60}")
        
        return final_M_matrices, final_A_matrices
    
    def _enforce_sign_convention(self, component: np.ndarray) -> np.ndarray:
        """Ensure component's largest absolute value is positive."""
        component = component.copy()
        max_idx = np.argmax(np.abs(component))
        if component[max_idx] < 0:
            component *= -1
        return component
    
    def _cluster_core_components(
        self, 
        core_components: Dict[str, List[np.ndarray]], 
        species_list: List[str], 
        num_runs: int,
        n_species: int
    ) -> Dict:
        """Cluster core components across all species."""
        results = {'clusters': [], 'centroids': {}}
        
        # Calculate scaled parameters
        base_per_dataset = 50
        base_runs = 100
        
        # Scale based on both num_runs and num_datasets
        runs_scaling_factor = num_runs / base_runs
        core_min_cluster_size = max(2, int(base_per_dataset * n_species * runs_scaling_factor))
        core_min_samples = max(2, int(base_per_dataset * n_species * runs_scaling_factor))
        core_min_per_dataset_count = max(1, int(base_per_dataset * runs_scaling_factor))
        
        # Prepare data matrix and dataset labels
        X_list, dataset_labels = [], []
        for dataset_idx, species in enumerate(species_list):
            if not core_components[species]:
                continue
            X = np.unique(np.array(core_components[species]), axis=0)
            X_list.append(X)
            dataset_labels.append(np.full(X.shape[0], dataset_idx))
        
        if not X_list:
            return results
        
        X = np.vstack(X_list)
        dataset_labels = np.concatenate(dataset_labels)
        
        # Cluster components
        clusterer = HDBSCAN(
            min_cluster_size=core_min_cluster_size,
            min_samples=core_min_samples,
            cluster_selection_epsilon=0.0,
            metric='euclidean',
            n_jobs=-1
        )
        labels = clusterer.fit_predict(X)
        
        # Identify valid clusters
        valid_clusters = []
        for label in np.unique(labels):
            if label == -1:
                continue
            
            mask = labels == label
            cluster_size = mask.sum()
            cluster_datasets = dataset_labels[mask]
            dataset_counts = np.bincount(cluster_datasets, minlength=len(species_list))
            
            # Check if cluster meets requirements
            if cluster_size >= core_min_cluster_size and np.all(dataset_counts >= core_min_per_dataset_count):
                valid_clusters.append({
                    'label': label,
                    'size': cluster_size,
                    'dataset_counts': dataset_counts
                })
        
        # Extract centroids
        for cluster_idx, cluster in enumerate(valid_clusters):
            cluster_label = cluster['label']
            cluster_centroids = {}
            
            for dataset_idx, species in enumerate(species_list):
                if cluster['dataset_counts'][dataset_idx] >= core_min_per_dataset_count:
                    dataset_mask = (labels == cluster_label) & (dataset_labels == dataset_idx)
                    centroid = np.mean(X[dataset_mask], axis=0)
                    cluster_centroids[species] = centroid
            
            if cluster_centroids:
                results['clusters'].append(cluster)
                results['centroids'][f"Core_{cluster_idx + 1}"] = cluster_centroids
        
        return results
    
    def _cluster_unique_components(
        self,
        unique_components: Dict[str, List[np.ndarray]],
        species_list: List[str],
        num_runs: int
    ) -> Dict[str, List[Tuple[str, np.ndarray]]]:
        """Cluster unique components per species."""
        results = {species: [] for species in species_list}
        
        # Calculate scaled parameters
        base_unique_min = 50
        base_unique_max = 100
        base_runs = 100
        
        # Scale based on num_runs only
        runs_scaling_factor = num_runs / base_runs
        unique_min_cluster_size = max(2, int(base_unique_min * runs_scaling_factor))
        unique_min_samples = max(2, int(base_unique_min * runs_scaling_factor))
        unique_max_cluster_size = int(base_unique_max * runs_scaling_factor)
        
        for species in species_list:
            components = unique_components[species]
            if not components:
                continue
            
            X = np.unique(np.array(components), axis=0)
            if X.size == 0:
                continue
            
            # Cluster components
            clusterer = HDBSCAN(
                min_cluster_size=unique_min_cluster_size,
                min_samples=unique_min_samples,
                cluster_selection_epsilon=0.0,
                metric='euclidean',
                n_jobs=-1
            )
            labels = clusterer.fit_predict(X)
            
            # Extract valid clusters
            valid_clusters = []
            for label in np.unique(labels):
                if label == -1:
                    continue
                mask = labels == label
                cluster_size = mask.sum()
                if unique_min_cluster_size <= cluster_size <= unique_max_cluster_size:
                    centroid = np.mean(X[mask], axis=0)
                    valid_clusters.append((f"Unique_{len(valid_clusters)+1}", centroid))
            
            results[species] = valid_clusters
        
        return results
    
    # Delegate methods to their respective modules
    def generate_BBH(self, output_path: str = "Output_BBH", threads: int = 1):
        """Generate BBH files using existing protein.faa files from each strain."""
        return generate_BBH(self, output_path, threads)
    
    def align_genes(self, input_bbh_dir: str = "Output_BBH", output_dir: str = "Output_Gene_Info", 
                    reference_order: Optional[List[str]] = None, bbh_threshold: Optional[float] = None) -> pd.DataFrame:
        """Align genes across all species using Union-Find algorithm."""
        return align_genes(self, input_bbh_dir, output_dir, reference_order, bbh_threshold)
    
    def optimize_number_of_core_components(self, **kwargs) -> Tuple[int, Dict[int, float]]:
        """Optimize the number of core components using specified metric."""
        return optimize_number_of_core_components(self, **kwargs)
    
    def optimize_number_of_unique_components(self, **kwargs) -> Tuple[Dict[str, int], Dict[str, int]]:
        """Optimize the number of unique components for each species."""
        return optimize_number_of_unique_components(self, **kwargs)
    
    def gff2pandas(self, gff_file: str, feature: str = "CDS", index: Optional[str] = None) -> pd.DataFrame:
        """Convert GFF file to pandas DataFrame."""
        return gff2pandas(gff_file, feature, index)
    
    def create_gene_table(self) -> None:
        """Create gene tables for all species."""
        return create_gene_table(self)
    
    def add_eggnog_annotation(self, eggnog_output_path: str) -> None:
        """
        Add eggNOG annotations to gene tables for all species.
        
        This method reads eggNOG-mapper output files and adds annotation columns
        to the existing gene tables. The annotations are matched using the 
        ncbi_protein column in gene_table.
        
        Parameters
        ----------
        eggnog_output_path : str
            Path to the eggNOG-mapper output directory containing subfolders
            for each species (e.g., ../imminer_2_industrial_strain/Output_eggnog_mapper)
            
        Notes
        -----
        The method expects the following structure:
        - eggnog_output_path/
            - species1/
                - species1.emapper.annotations
            - species2/
                - species2.emapper.annotations
            ...
            
        The following annotation columns will be added to gene_table:
        'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl',
        'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC', 
        'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction',
        'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'
        """
        import os
        
        eggnog_path = Path(eggnog_output_path)
        if not eggnog_path.exists():
            raise ValueError(f"eggNOG output path not found: {eggnog_output_path}")
        
        print(f"\nAdding eggNOG annotations from {eggnog_output_path}")
        print("=" * 60)
        
        # Columns to add from eggNOG annotation
        eggnog_columns = [
            'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl',
            'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC',
            'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction',
            'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'
        ]
        
        for species_name, species_data in self._species_data.items():
            print(f"\nProcessing {species_name}...")
            
            # Check if gene_table exists
            if species_data._gene_table is None:
                logger.warning(f"Gene table not found for {species_name}. Run create_gene_table() first.")
                print(f"  ✗ Gene table not found. Skipping...")
                continue
            
            # Find eggNOG annotation file
            species_eggnog_dir = eggnog_path / species_name
            if not species_eggnog_dir.exists():
                logger.warning(f"No eggNOG output directory found for {species_name}")
                print(f"  ✗ No eggNOG output directory found. Skipping...")
                continue
            
            annotation_file = species_eggnog_dir / f"{species_name}.emapper.annotations"
            if not annotation_file.exists():
                logger.warning(f"No annotation file found for {species_name}: {annotation_file}")
                print(f"  ✗ Annotation file not found. Skipping...")
                continue
            
            try:
                # Read eggNOG annotation file
                print(f"  - Reading {annotation_file.name}")
                df_eggnog = pd.read_csv(
                    annotation_file,
                    sep='\t',
                    skiprows=4,
                    dtype=str  # Read all columns as strings initially
                )
                
                # Remove the last 3 rows (summary lines)
                df_eggnog = df_eggnog[:-3]
                
                # Convert numeric columns
                df_eggnog['evalue'] = pd.to_numeric(df_eggnog['evalue'], errors='coerce')
                df_eggnog['score'] = pd.to_numeric(df_eggnog['score'], errors='coerce')
                
                # Rename query column to match ncbi_protein
                df_eggnog = df_eggnog.rename(columns={'#query': 'ncbi_protein'})
                
                # Get current gene table
                gene_table = species_data._gene_table.copy()
                
                # Check if ncbi_protein column exists in gene_table
                if 'ncbi_protein' not in gene_table.columns:
                    logger.warning(f"ncbi_protein column not found in gene_table for {species_name}")
                    print(f"  ✗ ncbi_protein column not found in gene_table. Skipping...")
                    continue
                
                # Reset index to have locus_tag as a column for merging
                gene_table_reset = gene_table.reset_index()
                
                # Merge with eggNOG annotations
                merged_table = gene_table_reset.merge(
                    df_eggnog[['ncbi_protein'] + eggnog_columns],
                    on='ncbi_protein',
                    how='left'
                )
                
                # Set locus_tag back as index
                merged_table.set_index('locus_tag', inplace=True)
                
                # Update species gene_table
                species_data._gene_table = merged_table
                
                # Count successful annotations
                annotated_count = merged_table[eggnog_columns[0]].notna().sum()
                total_genes = len(merged_table)
                
                print(f"  ✓ Added eggNOG annotations to {annotated_count}/{total_genes} genes")
                logger.info(f"eggNOG annotations added for {species_name}: {annotated_count}/{total_genes} genes")
                
            except Exception as e:
                logger.error(f"Error processing eggNOG annotations for {species_name}: {e}")
                print(f"  ✗ Error: {e}")
        
        print("\n" + "=" * 60)
        print("eggNOG annotation addition completed!")
    
    
    def optimize_M_thresholds(self, method: str = "Otsu's method", quantile_threshold: float = 90):
        """
        Optimize thresholds for M matrices across all species.
        
        This method calculates optimal thresholds for each component in the M matrices
        using the specified method, then creates M_thresholds and presence_matrix
        for each species.
        
        Parameters
        ----------
        method : str, optional
            Method to use for threshold optimization. Default is "Otsu's method".
            Currently only "Otsu's method" is supported.
        quantile_threshold : float, optional
            Percentile threshold for pre-filtering in Otsu's method. Default is 90.
            This removes the bottom X% of absolute values before applying Otsu's method
            to better identify true outliers in heavy-tailed distributions.
            
        Notes
        -----
        The Otsu's method implementation uses quantile-based pre-filtering to handle
        heavy-tailed distributions typical in gene expression data.
        """
        if method != "Otsu's method":
            raise ValueError(f"Method '{method}' not supported. Only 'Otsu's method' is currently available.")
        
        logger.info("Optimizing M thresholds using Otsu's method")
        
        # Check if M matrices exist
        for species_name, species_data in self._species_data.items():
            if species_data._M is None:
                raise ValueError(f"M matrix not found for {species_name}. Please run ICA first.")
        
        # Process each species
        for species_name, species_data in self._species_data.items():
            print(f"\nOptimizing thresholds for {species_name}...")
            
            M = species_data.M
            
            # Create M_thresholds dataframe
            thresholds = {}
            
            # Calculate threshold for each component
            for component in M.columns:
                threshold = self._quantile_otsu_threshold(M, component, quantile_threshold)
                thresholds[component] = threshold
            
            # Create M_thresholds dataframe with component names as index
            M_thresholds_df = pd.DataFrame(
                thresholds.values(),
                index=thresholds.keys(),
                columns=['M_threshold']
            )
            
            # Create presence matrix
            presence_matrix = pd.DataFrame(
                index=M.index,
                columns=M.columns,
                dtype=int
            )
            
            # Binarize based on thresholds
            for component in M.columns:
                threshold = thresholds[component]
                # Set to 1 if absolute value is above threshold, 0 otherwise
                presence_matrix[component] = (M[component].abs() > threshold).astype(int)
            
            # Store results
            species_data._M_thresholds = M_thresholds_df
            species_data._presence_matrix = presence_matrix
            
            # Report statistics
            n_components = len(M.columns)
            avg_genes_per_component = presence_matrix.sum().mean()
            print(f"✓ Optimized thresholds for {n_components} components")
            print(f"  Average genes per component: {avg_genes_per_component:.1f}")
        
        print("\nThreshold optimization completed!")
    
    def _quantile_otsu_threshold(self, multiModulon_M: pd.DataFrame, component_name: str, quantile_threshold: float = 90) -> float:
        """
        Calculate Otsu threshold using quantile-based pre-filtering.
        This method removes the central mass of data before applying Otsu,
        which helps identify true outliers in heavy-tailed distributions.
        
        Parameters
        ----------
        multiModulon_M : pd.DataFrame
            DataFrame with gene weights (genes as rows, components as columns)
        component_name : str
            Name of the component column to analyze
        quantile_threshold : float, optional
            Percentile threshold for pre-filtering. Default is 90.
        
        Returns
        -------
        threshold : float
            Optimal threshold value for the component
        """
        import numpy as np
        
        # Get weights for the specified component
        weights = multiModulon_M[component_name].values
        abs_weights = np.abs(weights)
        
        # Pre-filter: remove the bottom X% of absolute values
        pre_percentile = quantile_threshold
        pre_threshold = np.percentile(abs_weights, pre_percentile)
        filtered_weights = abs_weights[abs_weights > pre_threshold]
        
        # If too few points remain, return the percentile threshold
        if len(filtered_weights) < 10:
            return pre_threshold
        
        # Apply Otsu to the filtered data
        n_bins = 500
        counts, bin_edges = np.histogram(filtered_weights, bins=n_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Otsu algorithm on filtered data
        best_threshold = pre_threshold
        max_variance = 0
        
        total = counts.sum()
        if total == 0:
            return pre_threshold
            
        sum_total = (counts * bin_centers).sum()
        
        w0 = 0
        sum0 = 0
        
        for i in range(len(counts)):
            w0 += counts[i]
            if w0 == 0:
                continue
                
            w1 = total - w0
            if w1 == 0:
                break
                
            sum0 += counts[i] * bin_centers[i]
            m0 = sum0 / w0
            m1 = (sum_total - sum0) / w1
            
            # Between-class variance
            variance = w0 * w1 * (m0 - m1) ** 2
            
            if variance > max_variance:
                max_variance = variance
                best_threshold = bin_centers[i]
        
        return best_threshold
    
    def calculate_explained_variance(self) -> Dict[str, float]:
        """
        Calculate the explained variance for each species using X, M, and A matrices.
        
        This method computes the explained variance based on the reconstruction
        of the expression matrix X from the product of M and A matrices.
        The algorithm follows the approach of calculating individual component
        contributions and summing their explained variances.
        
        Returns
        -------
        dict
            Dictionary mapping species names to their total explained variance values.
            Keys are species names, values are the sum of explained variances.
            
        Raises
        ------
        ValueError
            If X, M, or A matrices are not available for any species.
            
        Examples
        --------
        >>> # After running multi-view ICA and generating A matrices
        >>> explained_var = multiModulon.calculate_explained_variance()
        >>> print(explained_var)
        {'strain1': 0.85, 'strain2': 0.82, 'strain3': 0.87}
        """
        explained_variance_results = {}
        
        # Check if ICA has been run
        ica_run = False
        for species_data in self._species_data.values():
            if species_data.M is not None:
                ica_run = True
                break
        
        if not ica_run:
            warnings.warn(
                "M and A matrices have not been generated from run_multiview_ica() or "
                "run_robust_multiview_ica(). Please run one of these methods first to "
                "ensure valid M and A matrices are available for explained variance calculation.",
                UserWarning
            )
        
        for species_name, species_data in self._species_data.items():
            # Check if all required matrices are available
            if species_data.X is None:
                raise ValueError(f"X matrix not found for {species_name}. Please run generate_X() first.")
            if species_data.M is None:
                raise ValueError(f"M matrix not found for {species_name}. Please run multi-view ICA first.")
            if species_data.A is None:
                raise ValueError(f"A matrix not found for {species_name}. Please run generate_A() first.")
            
            # Get matrices
            X = species_data.X
            M = species_data.M
            A = species_data.A
            
            # Center the data (following the algorithm)
            centered = X
            baseline = centered.subtract(centered.mean(axis=0), axis=1)
            
            # Initialize variables
            base_err = np.linalg.norm(baseline) ** 2
            MA = np.zeros(baseline.shape)
            rec_var = [0]
            ma_arrs = {}
            ma_weights = {}
            explained_variance_dict = {}
            
            # Get genes, samples, and components
            genes = X.index
            samples = X.columns
            imodulons = M.columns
            
            # Calculate individual component contributions
            for k in imodulons:
                # Compute outer product for this component
                ma_arr = np.dot(
                    M.loc[genes, k].values.reshape(len(genes), 1),
                    A.loc[k, samples].values.reshape(1, len(samples)),
                )
                ma_arrs[k] = ma_arr
                ma_weights[k] = np.sum(ma_arr**2)
            
            # Sort components by importance (highest weight first)
            sorted_mods = sorted(ma_weights, key=ma_weights.get, reverse=True)
            
            # Compute reconstructed variance
            i = 0
            for k in sorted_mods:
                MA = MA + ma_arrs[k]
                sa_err = np.linalg.norm(MA - baseline) ** 2
                rec_var.append((1 - sa_err / base_err))
                explained_variance_dict[k] = rec_var[i+1] - rec_var[i]
                i += 1
            
            # Sum all explained variances for this species
            total_explained_variance = sum(explained_variance_dict.values())
            explained_variance_results[species_name] = total_explained_variance
            
            # Log the result
            logger.info(f"Explained variance for {species_name}: {total_explained_variance:.4f}")
        
        return explained_variance_results
    
    def view_iModulon_weights(self, species: str, component: str, save_path: Optional[str] = None, 
                      fig_size: Tuple[float, float] = (6, 4), font_path: Optional[str] = None,
                      show_COG: bool = False):
        """
        Visualize gene weights for a specific iModulon component in a species.
        
        Creates a scatter plot showing gene weights across the genome, with genes
        positioned by their genomic coordinates.
        
        Parameters
        ----------
        species : str
            Species/strain name
        component : str
            Component name (e.g., 'Core_1', 'Unique_1')
        save_path : str, optional
            Path to save the plot. Can be a directory or file path.
            If directory, saves as '{species}_{component}_iModulon.svg'
            If None, displays plot without saving
        fig_size : tuple, optional
            Figure size as (width, height). Default: (6, 4)
        font_path : str, optional
            Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
            If provided, uses this font for all text elements
        show_COG : bool, optional
            If True, color genes above threshold by their COG category. Default: False
            
        Raises
        ------
        ValueError
            If species not found, M matrix not available, or component not found
        """
        # Validate species
        if species not in self._species_data:
            raise ValueError(f"Species '{species}' not found in loaded data")
        
        species_data = self._species_data[species]
        
        # Check if M matrix exists
        if species_data.M is None:
            raise ValueError(f"M matrix not found for {species}. Please run ICA first.")
        
        # Check if component exists
        if component not in species_data.M.columns:
            raise ValueError(f"Component '{component}' not found in M matrix. "
                           f"Available components: {list(species_data.M.columns)}")
        
        # Check if gene_table exists
        if species_data.gene_table is None:
            raise ValueError(f"Gene table not found for {species}. Please run create_gene_table() first.")
        
        # Get gene weights for the component
        gene_weights = species_data.M[component]
        
        # Get gene positions from gene_table
        gene_table = species_data.gene_table
        
        # Create mapping from leftmost genes (in M matrix) to species genes (in gene_table)
        leftmost_to_species = {}
        species_to_leftmost = {}
        
        if self.combined_gene_db is not None and species in self.combined_gene_db.columns:
            # First, find the leftmost gene for each row (same logic as in alignment)
            for idx, row in self.combined_gene_db.iterrows():
                # Find the leftmost (first non-None) gene in this row
                leftmost_gene = None
                for col in self.combined_gene_db.columns:
                    val = row[col]
                    if pd.notna(val) and val != "None" and val is not None:
                        leftmost_gene = val
                        break
                
                # Get the species-specific gene
                species_gene = row[species]
                
                if leftmost_gene and pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                    leftmost_to_species[leftmost_gene] = species_gene
                    species_to_leftmost[species_gene] = leftmost_gene
        
        # Create lists for plotting
        genes_for_plotting = []
        weights_for_plotting = []
        positions_for_plotting = []
        
        # Debug: log mapping statistics
        logger.debug(f"Total genes in M matrix: {len(gene_weights)}")
        logger.debug(f"Leftmost to species mappings found: {len(leftmost_to_species)}")
        
        # Process each gene in M matrix (which uses leftmost names)
        genes_not_found = []
        for gene_in_M in gene_weights.index:
            weight = gene_weights[gene_in_M]
            
            # Find the corresponding species gene name
            if gene_in_M in leftmost_to_species:
                # This is a leftmost gene, map to species gene
                species_gene = leftmost_to_species[gene_in_M]
                if species_gene in gene_table.index:
                    gene_info = gene_table.loc[species_gene]
                    center = (gene_info['start'] + gene_info['end']) / 2
                    genes_for_plotting.append(species_gene)
                    weights_for_plotting.append(weight)
                    positions_for_plotting.append(center / 1e6)  # Convert to Mb
                else:
                    genes_not_found.append(f"{gene_in_M}->{species_gene} (not in gene_table)")
            elif gene_in_M in gene_table.index:
                # This gene name exists directly in gene_table (not renamed)
                gene_info = gene_table.loc[gene_in_M]
                center = (gene_info['start'] + gene_info['end']) / 2
                genes_for_plotting.append(gene_in_M)
                weights_for_plotting.append(weight)
                positions_for_plotting.append(center / 1e6)  # Convert to Mb
            else:
                genes_not_found.append(f"{gene_in_M} (no mapping found)")
        
        if genes_not_found:
            logger.debug(f"Genes not plotted ({len(genes_not_found)}): {genes_not_found[:10]}...")
            logger.debug(f"Genes successfully plotted: {len(genes_for_plotting)}/{len(gene_weights)}")
        
        if not genes_for_plotting:
            raise ValueError(f"No gene position information found for genes in {species}")
        
        # Now we have the correct mapping
        x_positions = positions_for_plotting
        y_weights = weights_for_plotting
        genes_with_pos = genes_for_plotting
        
        # Get threshold if available
        threshold = None
        if hasattr(species_data, '_M_thresholds') and species_data._M_thresholds is not None:
            if component in species_data._M_thresholds.index:
                threshold = species_data._M_thresholds.loc[component, 'M_threshold']
        
        # Create figure - the given fig_size is for the plot area only
        # If showing COG, we need extra space for the legend
        if show_COG and 'COG_category' in gene_table.columns:
            # Add extra width for legend (80% more for long COG category names)
            fig_width = fig_size[0] * 1.8
            fig = plt.figure(figsize=(fig_width, fig_size[1]))
            # Create axis that uses only the original fig_size width
            ax = fig.add_axes([0.1, 0.1, fig_size[0]/fig_width * 0.8, 0.8])
        else:
            fig, ax = plt.subplots(figsize=fig_size)
        
        # Set font properties if provided
        if font_path and os.path.exists(font_path):
            font_prop = fm.FontProperties(fname=font_path)
            plt.rcParams['font.family'] = font_prop.get_name()
        
        # Create scatter plot with color coding
        if show_COG and 'COG_category' in gene_table.columns and 'Description' in gene_table.columns:
            # Helper function to get COG color
            def get_cog_color(gene_name, weight):
                if threshold is not None and abs(weight) <= threshold:
                    return 'grey'  # Below threshold
                
                # Try to get COG category for this gene
                cog_cat = None
                if gene_name in gene_table.index:
                    cog_info = gene_table.loc[gene_name, 'COG_category']
                    desc_info = gene_table.loc[gene_name, 'Description']
                    
                    # Check COG_category first
                    if pd.notna(cog_info) and cog_info != '-':
                        # Extract first letter if it's a letter code
                        if isinstance(cog_info, str) and len(cog_info) > 0:
                            first_letter = cog_info[0]
                            if first_letter in COG_LETTER_CODES:
                                cog_cat = COG_LETTER_CODES[first_letter]
                    
                    # If not found in COG_category, check Description
                    if cog_cat is None and pd.notna(desc_info) and desc_info != '-':
                        # Check if description matches any COG category
                        for cat_name in COG_COLORS.keys():
                            if cat_name != 'No COG annotation' and cat_name.lower() in str(desc_info).lower():
                                cog_cat = cat_name
                                break
                
                if cog_cat and cog_cat in COG_COLORS:
                    return COG_COLORS[cog_cat]
                else:
                    return COG_COLORS['No COG annotation']
            
            # Get colors for all genes
            colors = [get_cog_color(gene, weight) for gene, weight in zip(genes_with_pos, y_weights)]
            
            # Create scatter plot
            scatter = ax.scatter(x_positions, y_weights, alpha=0.6, s=20, c=colors)
            
            # Add horizontal threshold lines if available
            if threshold is not None:
                ax.axhline(y=threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
                ax.axhline(y=-threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
            
            # Create legend
            unique_colors = {}
            for gene, color, weight in zip(genes_with_pos, colors, y_weights):
                if threshold is None or abs(weight) > threshold:
                    # Get COG category name for this color
                    for cat_name, cat_color in COG_COLORS.items():
                        if cat_color == color:
                            unique_colors[cat_name] = color
                            break
            
            # Sort categories by type
            info_storage = ['Translation, ribosomal structure and biogenesis', 'RNA processing and modification', 
                           'Transcription', 'Replication, recombination and repair', 'Chromatin structure and dynamics']
            cellular = ['Cell cycle control, cell division, chromosome partitioning', 'Nuclear structure', 
                       'Defense mechanisms', 'Signal transduction mechanisms', 'Cell wall/membrane/envelope biogenesis',
                       'Cell motility', 'Cytoskeleton', 'Extracellular structures', 
                       'Intracellular trafficking, secretion, and vesicular transport',
                       'Posttranslational modification, protein turnover, chaperones',
                       'Post-translational modification, protein turnover, and chaperones']
            metabolism = ['Energy production and conversion', 'Carbohydrate transport and metabolism',
                         'Amino acid transport and metabolism', 'Nucleotide transport and metabolism',
                         'Coenzyme transport and metabolism', 'Lipid transport and metabolism',
                         'Inorganic ion transport and metabolism', 
                         'Secondary metabolites biosynthesis, transport and catabolism',
                         'Secondary metabolites biosynthesis, transport, and catabolism']
            poorly_char = ['Function unknown', 'No COG annotation']
            
            sorted_categories = []
            for cat_list in [info_storage, cellular, metabolism, poorly_char]:
                for cat in cat_list:
                    if cat in unique_colors:
                        sorted_categories.append(cat)
            
            if sorted_categories:
                # Create legend elements
                from matplotlib.patches import Patch
                legend_elements = [Patch(facecolor=unique_colors[cat], label=cat) 
                                 for cat in sorted_categories]
                
                # Add legend with appropriate font
                if font_path and os.path.exists(font_path):
                    legend = ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5),
                                     frameon=True, fontsize=10)
                    for text in legend.get_texts():
                        text.set_fontproperties(font_prop)
                else:
                    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5),
                            frameon=True, fontsize=10)
        
        elif threshold is not None:
            # Original threshold-based coloring
            colors = []
            for weight in y_weights:
                if abs(weight) > threshold:
                    colors.append('lightblue')  # Above threshold
                else:
                    colors.append('grey')  # Below threshold
            ax.scatter(x_positions, y_weights, alpha=0.6, s=20, c=colors)
            
            # Add horizontal threshold lines
            ax.axhline(y=threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
            ax.axhline(y=-threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
        else:
            # No threshold available, use default plotting
            ax.scatter(x_positions, y_weights, alpha=0.6, s=20)
        
        # Set labels and title
        ax.set_xlabel('Gene Start (1e6)', fontsize=12)
        ax.set_ylabel('Gene Weight', fontsize=12)
        ax.set_title(f'iModulon {component} on {species}', fontsize=14)
        
        # Set font for tick labels if font_path provided
        if font_path and os.path.exists(font_path):
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontproperties(font_prop)
            ax.xaxis.label.set_fontproperties(font_prop)
            ax.yaxis.label.set_fontproperties(font_prop)
            ax.title.set_fontproperties(font_prop)
        
        # Tight layout - only apply if not showing COG (manual layout for COG)
        if not (show_COG and 'COG_category' in gene_table.columns):
            plt.tight_layout()
        
        # Save or show
        if save_path:
            # Determine save file path
            save_path = Path(save_path)
            if save_path.suffix in ['.svg', '.png', '.pdf', '.jpg']:
                # Full file path provided
                save_file = save_path
            else:
                # Directory provided, use default name
                save_path.mkdir(parents=True, exist_ok=True)
                save_file = save_path / f"{species}_{component}_iModulon.svg"
            
            # Save with tight bbox to include legend if present
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_file}")
            plt.close()
        else:
            plt.show()
    
    def view_iModulon_activities(self, species: str, component: str, save_path: Optional[str] = None,
                                fig_size: Tuple[float, float] = (12, 3), font_path: Optional[str] = None,
                                highlight_project: Optional[Union[str, List[str]]] = None, highlight_study: Optional[str] = None):
        """
        Visualize iModulon activities for a specific component in a species.
        
        Creates a bar plot showing component activities across samples, with samples
        grouped by project/study_accession from the sample sheet.
        
        Parameters
        ----------
        species : str
            Species/strain name
        component : str
            Component name (e.g., 'Core_1', 'Unique_1')
        save_path : str, optional
            Path to save the plot. Can be a directory or file path.
            If directory, saves as '{species}_{component}_activities.svg'
            If None, displays plot without saving
        fig_size : tuple, optional
            Figure size as (width, height). Default: (12, 3)
        font_path : str, optional
            Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
            If provided, uses this font for all text elements
        highlight_project : str or list of str, optional
            Project name(s) to highlight with different colors. Uses 'condition' column for legend.
            Can be a single project name (string) or multiple project names (list)
        highlight_study : str, optional
            Study accession to highlight with different colors. Uses 'sample_description' column for legend
            
        Raises
        ------
        ValueError
            If species not found, A matrix not available, or component not found
        """
        # Validate species
        if species not in self._species_data:
            raise ValueError(f"Species '{species}' not found in loaded data")
        
        species_data = self._species_data[species]
        
        # Check if A matrix exists
        if species_data.A is None:
            raise ValueError(f"A matrix not found for {species}. Please run ICA first.")
        
        # Check if component exists
        if component not in species_data.A.index:
            raise ValueError(f"Component '{component}' not found in A matrix. "
                           f"Available components: {list(species_data.A.index)}")
        
        # Get component activities
        activities = species_data.A.loc[component]
        
        # Get sample sheet
        sample_sheet = species_data.sample_sheet
        condition_mode = False
        condition_data = {}
        
        if sample_sheet is None:
            # If no sample sheet, use simple bar plot without grouping
            x_positions = range(len(activities))
            x_labels = activities.index.tolist()
            group_labels = []
        else:
            # Check if condition column exists
            if 'condition' in sample_sheet.columns:
                condition_mode = True
                # Ensure activities and sample_sheet are aligned
                common_samples = [s for s in activities.index if s in sample_sheet.index]
                activities = activities[common_samples]
                sample_sheet = sample_sheet.loc[common_samples]
                
                # Group samples by condition
                for sample in common_samples:
                    condition = sample_sheet.loc[sample, 'condition']
                    if condition not in condition_data:
                        condition_data[condition] = {
                            'samples': [],
                            'activities': [],
                            'mean_activity': 0
                        }
                    condition_data[condition]['samples'].append(sample)
                    condition_data[condition]['activities'].append(activities[sample])
                
                # Calculate mean activities for each condition
                conditions = list(condition_data.keys())
                mean_activities = []
                for condition in conditions:
                    mean_act = np.mean(condition_data[condition]['activities'])
                    condition_data[condition]['mean_activity'] = mean_act
                    mean_activities.append(mean_act)
                
                x_positions = range(len(conditions))
                x_labels = conditions
                group_labels = []
            else:
                # Original logic for project/study grouping
                # Determine grouping column
                group_col = None
                if 'project' in sample_sheet.columns:
                    group_col = 'project'
                elif 'study_accession' in sample_sheet.columns:
                    group_col = 'study_accession'
                
                if group_col:
                    # Ensure activities and sample_sheet are aligned
                    common_samples = [s for s in activities.index if s in sample_sheet.index]
                    activities = activities[common_samples]
                    sample_sheet = sample_sheet.loc[common_samples]
                    
                    # Get group labels
                    group_labels = sample_sheet[group_col].fillna('Unknown').astype(str).tolist()
                    unique_groups = []
                    for g in group_labels:
                        if g not in unique_groups:
                            unique_groups.append(g)
                    
                    x_positions = range(len(activities))
                    x_labels = [unique_groups.index(g) if g in unique_groups else -1 for g in group_labels]
                else:
                    # No grouping column available
                    x_positions = range(len(activities))
                    x_labels = activities.index.tolist()
                    group_labels = []
        
        # Prepare colors for bars
        colors = []
        legend_elements = []
        
        if condition_mode:
            # In condition mode, use default colors unless highlighting is specified
            colors = ['lightblue'] * len(conditions)
        elif highlight_project and sample_sheet is not None and 'project' in sample_sheet.columns:
            # Convert single string to list for uniform handling
            if isinstance(highlight_project, str):
                highlight_projects = [highlight_project]
            else:
                highlight_projects = highlight_project
            
            # Color by project with condition as legend
            for i, sample in enumerate(activities.index):
                if sample in sample_sheet.index:
                    project = str(sample_sheet.loc[sample, 'project'])
                    # Check if project is in the highlight list
                    project_highlighted = any(str(hp) == project for hp in highlight_projects)
                    
                    if project_highlighted:
                        condition = sample_sheet.loc[sample, 'condition'] if 'condition' in sample_sheet.columns else 'Unknown'
                        # Assign unique color per condition
                        if condition not in [elem[0] for elem in legend_elements]:
                            color = plt.cm.tab10(len(legend_elements))
                            legend_elements.append((condition, color))
                            colors.append(color)
                        else:
                            # Find existing color for this condition
                            for cond, col in legend_elements:
                                if cond == condition:
                                    colors.append(col)
                                    break
                    else:
                        colors.append('lightblue')
                else:
                    colors.append('lightblue')
        elif highlight_study and sample_sheet is not None and 'study_accession' in sample_sheet.columns:
            # Color by study with sample_description as legend
            for i, sample in enumerate(activities.index):
                if sample in sample_sheet.index:
                    study = sample_sheet.loc[sample, 'study_accession']
                    if str(study) == str(highlight_study):
                        desc = sample_sheet.loc[sample, 'sample_description'] if 'sample_description' in sample_sheet.columns else 'Unknown'
                        # Assign unique color per description
                        if desc not in [elem[0] for elem in legend_elements]:
                            color = plt.cm.tab10(len(legend_elements))
                            legend_elements.append((desc, color))
                            colors.append(color)
                        else:
                            # Find existing color for this description
                            for d, col in legend_elements:
                                if d == desc:
                                    colors.append(col)
                                    break
                    else:
                        colors.append('lightblue')
                else:
                    colors.append('lightblue')
        else:
            # Default light blue color for all bars
            colors = ['lightblue'] * len(activities)
        
        # Create figure
        fig, ax = plt.subplots(figsize=fig_size)
        
        # Set font properties if provided
        if font_path and os.path.exists(font_path):
            font_prop = fm.FontProperties(fname=font_path)
            plt.rcParams['font.family'] = font_prop.get_name()
        
        # Create bar plot
        if condition_mode:
            # Plot averaged bars for conditions
            bars = ax.bar(x_positions, mean_activities, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
            
            # Add individual sample points as black dots
            for i, condition in enumerate(conditions):
                sample_activities = condition_data[condition]['activities']
                # Create small random jitter for x-position to avoid overlapping dots
                jitter = np.random.uniform(-0.2, 0.2, size=len(sample_activities))
                x_points = [i + j for j in jitter]
                ax.scatter(x_points, sample_activities, color='black', s=30, zorder=10, alpha=0.7)
        else:
            # Original bar plot for individual samples
            bars = ax.bar(x_positions, activities.values, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # Add horizontal line at y=0
        ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
        
        # Remove grid
        ax.grid(False)
        
        # Add vertical dotted lines to separate projects
        if group_col and group_labels:
            current_group = group_labels[0]
            for i in range(1, len(group_labels)):
                if group_labels[i] != current_group:
                    # Add vertical line at boundary
                    ax.axvline(x=i-0.5, color='gray', linestyle=':', linewidth=1, alpha=0.7)
                    current_group = group_labels[i]
        
        # Set labels and title
        if condition_mode:
            ax.set_xlabel('Conditions', fontsize=12)
        else:
            ax.set_xlabel('Samples', fontsize=12)
        ax.set_ylabel('iModulon Activity', fontsize=12)
        ax.set_title(f'iModulon {component} on {species}', fontsize=14)
        
        # Set x-axis ticks
        if condition_mode:
            # Show condition labels on x-axis
            ax.set_xticks(x_positions)
            ax.set_xticklabels(x_labels, rotation=45, ha='right')
        elif group_col and group_labels:
            # Show group labels on x-axis
            tick_positions = []
            tick_labels = []
            current_group = None
            group_start = 0
            
            for i, group in enumerate(group_labels):
                if group != current_group:
                    if current_group is not None:
                        # Add tick for previous group
                        tick_positions.append((group_start + i - 1) / 2)
                        tick_labels.append(current_group)
                    current_group = group
                    group_start = i
            
            # Add last group
            if current_group is not None:
                tick_positions.append((group_start + len(group_labels) - 1) / 2)
                tick_labels.append(current_group)
            
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, rotation=45, ha='right')
        else:
            # No grouping - don't show individual sample names
            ax.set_xticks([])
        
        # Add legend if highlighting
        if legend_elements:
            from matplotlib.patches import Patch
            patches = [Patch(color=color, label=label) for label, color in legend_elements]
            ax.legend(handles=patches, loc='center left', bbox_to_anchor=(1, 0.5))
            
            # Apply font to legend if provided
            if font_path and os.path.exists(font_path):
                for text in ax.get_legend().get_texts():
                    text.set_fontproperties(font_prop)
        
        # Set font for tick labels if font_path provided
        if font_path and os.path.exists(font_path):
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontproperties(font_prop)
            ax.xaxis.label.set_fontproperties(font_prop)
            ax.yaxis.label.set_fontproperties(font_prop)
            ax.title.set_fontproperties(font_prop)
        
        # Tight layout - adjust if legend is present
        if legend_elements:
            plt.tight_layout()
            # Make room for legend on the right
            plt.subplots_adjust(right=0.85)
        else:
            plt.tight_layout()
        
        # Save or show
        if save_path:
            # Determine save file path
            save_path = Path(save_path)
            if save_path.suffix in ['.svg', '.png', '.pdf', '.jpg']:
                # Full file path provided
                save_file = save_path
            else:
                # Directory provided, use default name
                save_path.mkdir(parents=True, exist_ok=True)
                save_file = save_path / f"{species}_{component}_activities.svg"
            
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_file}")
            plt.close()
        else:
            plt.show()
    
    def view_iModulon_genes(self, species: str, component: str) -> pd.DataFrame:
        """
        Return a subset of the gene table for genes in a specific iModulon component.
        
        Uses the presence_matrix to identify which genes belong to the specified component,
        then returns the corresponding rows from the gene_table.
        
        Parameters
        ----------
        species : str
            Species/strain name
        component : str
            Component name (e.g., 'Core_1', 'Unique_1')
            
        Returns
        -------
        pd.DataFrame
            Subset of gene_table containing only genes in the specified component.
            Returns empty DataFrame if no genes found or if data not available.
            
        Raises
        ------
        ValueError
            If species not found or if presence_matrix/gene_table not available
            
        Examples
        --------
        >>> gene_subset = mm.view_iModulon_genes('E_coli', 'Core_1')
        >>> print(f"Found {len(gene_subset)} genes in Core_1")
        """
        # Validate species
        if species not in self._species_data:
            raise ValueError(f"Species '{species}' not found in loaded data")
        
        species_data = self._species_data[species]
        
        # Check if presence_matrix exists
        if species_data._presence_matrix is None:
            raise ValueError(f"presence_matrix not found for {species}. Please run optimize_M_thresholds() first.")
        
        # Check if gene_table exists
        if species_data.gene_table is None:
            raise ValueError(f"gene_table not found for {species}.")
        
        # Check if component exists
        if component not in species_data._presence_matrix.columns:
            raise ValueError(f"Component '{component}' not found in presence_matrix. "
                           f"Available components: {list(species_data._presence_matrix.columns)}")
        
        # Get genes in this component from presence_matrix
        component_genes = species_data._presence_matrix[
            species_data._presence_matrix[component] == 1
        ].index.tolist()
        
        if not component_genes:
            logger.warning(f"No genes found in component '{component}' for species '{species}'")
            return pd.DataFrame()
        
        # Map component genes to gene_table
        # First check if we need to map through combined_gene_db
        genes_in_table = []
        
        if self.combined_gene_db is not None and species in self.combined_gene_db.columns:
            # Map from leftmost genes to species-specific genes
            for gene in component_genes:
                # Check if gene exists directly in gene_table first
                if gene in species_data.gene_table.index:
                    genes_in_table.append(gene)
                else:
                    # Try to map through combined_gene_db
                    for idx, row in self.combined_gene_db.iterrows():
                        # Find leftmost gene
                        leftmost_gene = None
                        for col in self.combined_gene_db.columns:
                            val = row[col]
                            if pd.notna(val) and val != "None" and val is not None:
                                leftmost_gene = val
                                break
                        
                        # If this is our gene, get the species-specific name
                        if leftmost_gene == gene:
                            species_gene = row[species]
                            if pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                                if species_gene in species_data.gene_table.index:
                                    genes_in_table.append(species_gene)
                                    break
        else:
            # No combined_gene_db, just check directly
            genes_in_table = [g for g in component_genes if g in species_data.gene_table.index]
        
        # Get subset of gene_table
        if genes_in_table:
            gene_subset = species_data.gene_table.loc[genes_in_table].copy()
            
            # Sort by start position if the column exists
            if 'start' in gene_subset.columns:
                gene_subset = gene_subset.sort_values('start')
                
            logger.info(f"Found {len(gene_subset)} genes in component '{component}' for species '{species}'")
            return gene_subset
        else:
            logger.warning(f"No genes from component '{component}' found in gene_table for species '{species}'")
            return pd.DataFrame()
    
    def view_core_iModulon_weights(self, component: str, save_path: Optional[str] = None,
                           fig_size: Tuple[float, float] = (6, 4), font_path: Optional[str] = None,
                           show_COG: bool = False, reference_order: Optional[List[str]] = None):
        """
        Visualize a core iModulon component across all species.
        
        Creates individual plots for each species showing the same core component,
        and saves them to the specified directory. This method calls view_iModulons
        for each species with the given core component.
        
        Parameters
        ----------
        component : str
            Core component name (e.g., 'Core_1', 'Core_2')
        save_path : str, optional
            Directory path to save all plots. If None, displays plots without saving.
            Each plot is saved as '{species}_{component}_iModulon.svg'
        fig_size : tuple, optional
            Figure size for individual plots as (width, height). Default: (6, 4)
        font_path : str, optional
            Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
            If provided, uses this font for all text elements in all plots
        show_COG : bool, optional
            If True, color genes above threshold by their COG category. Default: False
        reference_order : list of str, optional
            Custom order for arranging species in subplots when show_COG is True.
            First 3 species will be placed in the first row, remaining in the second row.
            Example: ['MG1655', 'BL21', 'C', 'Crooks', 'W', 'W3110']
            If not provided, species are plotted in their default order.
            
        Raises
        ------
        ValueError
            If component is not a core component or not found in any species
        """
        # Get all species
        species_list = list(self._species_data.keys())
        
        # Verify this is a core component by checking it exists in multiple species
        species_with_component = []
        for species in species_list:
            species_data = self._species_data[species]
            if species_data.M is not None and component in species_data.M.columns:
                species_with_component.append(species)
        
        if not species_with_component:
            raise ValueError(f"Component '{component}' not found in any species")
        
        if not component.startswith('Core_'):
            warnings.warn(f"Component '{component}' does not appear to be a core component "
                         f"(core components typically start with 'Core_')")
        
        # Create save directory if needed
        if save_path:
            save_path = Path(save_path)
            if not save_path.suffix:  # It's a directory
                save_path.mkdir(parents=True, exist_ok=True)
        
        # Generate plots for each species
        logger.info(f"Generating plots for core component '{component}' across {len(species_with_component)} species")
        
        if show_COG:
            # When showing COG, create a combined plot with subplots and single legend
            
            # Apply reference_order if provided
            if reference_order:
                # Filter reference_order to only include species that have the component
                ordered_species = [sp for sp in reference_order if sp in species_with_component]
                # Add any remaining species not in reference_order
                remaining_species = [sp for sp in species_with_component if sp not in ordered_species]
                species_with_component = ordered_species + remaining_species
            
            n_species = len(species_with_component)
            
            # Determine subplot layout
            if reference_order and n_species > 3:
                # Custom layout: first 3 in first row, rest in second row
                n_rows = 2
                n_cols = max(3, n_species - 3)
            elif n_species <= 3:
                n_rows = 1
                n_cols = n_species
            elif n_species <= 6:
                n_rows = 2
                n_cols = 3
            else:
                n_cols = 3
                n_rows = (n_species + 2) // 3
            
            # Create figure with subplots
            fig_width = fig_size[0] * n_cols
            fig_height = fig_size[1] * n_rows + 2  # Extra space for legend
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
            
            # Ensure axes is always a 2D array
            if n_species == 1:
                axes = np.array([[axes]])
            elif n_rows == 1:
                axes = axes.reshape(1, -1)
            elif n_cols == 1:
                axes = axes.reshape(-1, 1)
            
            # Set font properties if provided
            if font_path and os.path.exists(font_path):
                font_prop = fm.FontProperties(fname=font_path)
                plt.rcParams['font.family'] = font_prop.get_name()
            
            # Track all unique COG categories across all species
            all_unique_colors = {}
            
            # Plot each species
            for idx, species in enumerate(species_with_component):
                if reference_order and n_species > 3:
                    # Custom layout for reference_order
                    if idx < 3:
                        # First 3 species in first row
                        row = 0
                        col = idx
                    else:
                        # Remaining species in second row
                        row = 1
                        col = idx - 3
                else:
                    # Default layout
                    row = idx // n_cols
                    col = idx % n_cols
                ax = axes[row, col]
                
                species_data = self._species_data[species]
                
                # Get gene weights and gene table
                gene_weights = species_data.M[component]
                gene_table = species_data.gene_table
                
                # Create mapping from leftmost genes (in M matrix) to species genes (in gene_table)
                leftmost_to_species = {}
                
                if self.combined_gene_db is not None and species in self.combined_gene_db.columns:
                    # First, find the leftmost gene for each row (same logic as in alignment)
                    for idx, row in self.combined_gene_db.iterrows():
                        # Find the leftmost (first non-None) gene in this row
                        leftmost_gene = None
                        for col in self.combined_gene_db.columns:
                            val = row[col]
                            if pd.notna(val) and val != "None" and val is not None:
                                leftmost_gene = val
                                break
                        
                        # Get the species-specific gene
                        species_gene = row[species]
                        
                        if leftmost_gene and pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                            leftmost_to_species[leftmost_gene] = species_gene
                
                # Create lists for plotting
                genes_for_plotting = []
                weights_for_plotting = []
                positions_for_plotting = []
                
                # Process each gene in M matrix (which uses leftmost names)
                for gene_in_M in gene_weights.index:
                    weight = gene_weights[gene_in_M]
                    
                    # Find the corresponding species gene name
                    if gene_in_M in leftmost_to_species:
                        # This is a leftmost gene, map to species gene
                        species_gene = leftmost_to_species[gene_in_M]
                        if species_gene in gene_table.index:
                            gene_info = gene_table.loc[species_gene]
                            center = (gene_info['start'] + gene_info['end']) / 2
                            genes_for_plotting.append(species_gene)
                            weights_for_plotting.append(weight)
                            positions_for_plotting.append(center / 1e6)  # Convert to Mb
                    elif gene_in_M in gene_table.index:
                        # This gene name exists directly in gene_table (not renamed)
                        gene_info = gene_table.loc[gene_in_M]
                        center = (gene_info['start'] + gene_info['end']) / 2
                        genes_for_plotting.append(gene_in_M)
                        weights_for_plotting.append(weight)
                        positions_for_plotting.append(center / 1e6)  # Convert to Mb
                
                # Now we have the correct mapping
                x_positions = positions_for_plotting
                y_weights = weights_for_plotting
                genes_with_pos = genes_for_plotting
                
                # Get threshold
                threshold = None
                if hasattr(species_data, '_M_thresholds') and species_data._M_thresholds is not None:
                    if component in species_data._M_thresholds.index:
                        threshold = species_data._M_thresholds.loc[component, 'M_threshold']
                
                # Helper function to get COG color
                def get_cog_color(gene_name, weight):
                    if threshold is not None and abs(weight) <= threshold:
                        return 'grey'
                    
                    cog_cat = None
                    if gene_name in gene_table.index and 'COG_category' in gene_table.columns:
                        cog_info = gene_table.loc[gene_name, 'COG_category']
                        desc_info = gene_table.loc[gene_name, 'Description'] if 'Description' in gene_table.columns else None
                        
                        if pd.notna(cog_info) and cog_info != '-':
                            if isinstance(cog_info, str) and len(cog_info) > 0:
                                first_letter = cog_info[0]
                                if first_letter in COG_LETTER_CODES:
                                    cog_cat = COG_LETTER_CODES[first_letter]
                        
                        if cog_cat is None and desc_info and pd.notna(desc_info) and desc_info != '-':
                            for cat_name in COG_COLORS.keys():
                                if cat_name != 'No COG annotation' and cat_name.lower() in str(desc_info).lower():
                                    cog_cat = cat_name
                                    break
                    
                    if cog_cat and cog_cat in COG_COLORS:
                        return COG_COLORS[cog_cat]
                    else:
                        return COG_COLORS['No COG annotation']
                
                # Get colors
                colors = [get_cog_color(gene, weight) for gene, weight in zip(genes_with_pos, y_weights)]
                
                # Track unique colors for legend
                for gene, color, weight in zip(genes_with_pos, colors, y_weights):
                    if threshold is None or abs(weight) > threshold:
                        for cat_name, cat_color in COG_COLORS.items():
                            if cat_color == color:
                                all_unique_colors[cat_name] = color
                                break
                
                # Create scatter plot
                ax.scatter(x_positions, y_weights, alpha=0.6, s=20, c=colors)
                
                # Add threshold lines
                if threshold is not None:
                    ax.axhline(y=threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
                    ax.axhline(y=-threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
                
                # Set labels and title
                ax.set_xlabel('Gene Start (1e6)', fontsize=10)
                ax.set_ylabel('Gene Weight', fontsize=10)
                ax.set_title(f'{species}', fontsize=12)
                
                # Set font for labels if provided
                if font_path and os.path.exists(font_path):
                    for label in ax.get_xticklabels() + ax.get_yticklabels():
                        label.set_fontproperties(font_prop)
                    ax.xaxis.label.set_fontproperties(font_prop)
                    ax.yaxis.label.set_fontproperties(font_prop)
                    ax.title.set_fontproperties(font_prop)
            
            # Hide empty subplots
            for idx in range(n_species, n_rows * n_cols):
                row = idx // n_cols
                col = idx % n_cols
                axes[row, col].axis('off')
            
            # Add main title
            fig.suptitle(f'Core iModulon {component}', fontsize=16, y=0.98)
            if font_path and os.path.exists(font_path):
                fig.suptitle(f'Core iModulon {component}', fontsize=16, y=0.98, fontproperties=font_prop)
            
            # Create legend at bottom
            if all_unique_colors:
                # Sort categories
                info_storage = ['Translation, ribosomal structure and biogenesis', 'RNA processing and modification', 
                               'Transcription', 'Replication, recombination and repair', 'Chromatin structure and dynamics']
                cellular = ['Cell cycle control, cell division, chromosome partitioning', 'Nuclear structure', 
                           'Defense mechanisms', 'Signal transduction mechanisms', 'Cell wall/membrane/envelope biogenesis',
                           'Cell motility', 'Cytoskeleton', 'Extracellular structures', 
                           'Intracellular trafficking, secretion, and vesicular transport',
                           'Posttranslational modification, protein turnover, chaperones',
                           'Post-translational modification, protein turnover, and chaperones']
                metabolism = ['Energy production and conversion', 'Carbohydrate transport and metabolism',
                             'Amino acid transport and metabolism', 'Nucleotide transport and metabolism',
                             'Coenzyme transport and metabolism', 'Lipid transport and metabolism',
                             'Inorganic ion transport and metabolism', 
                             'Secondary metabolites biosynthesis, transport and catabolism',
                             'Secondary metabolites biosynthesis, transport, and catabolism']
                poorly_char = ['Function unknown', 'No COG annotation']
                
                sorted_categories = []
                for cat_list in [info_storage, cellular, metabolism, poorly_char]:
                    for cat in cat_list:
                        if cat in all_unique_colors:
                            sorted_categories.append(cat)
                
                # Create legend elements
                from matplotlib.patches import Patch
                legend_elements = [Patch(facecolor=all_unique_colors[cat], label=cat) 
                                 for cat in sorted_categories]
                
                # Add legend at bottom
                if font_path and os.path.exists(font_path):
                    legend = fig.legend(handles=legend_elements, loc='lower center', 
                                      bbox_to_anchor=(0.5, -0.03), ncol=3, frameon=True, fontsize=10)
                    for text in legend.get_texts():
                        text.set_fontproperties(font_prop)
                else:
                    fig.legend(handles=legend_elements, loc='lower center', 
                             bbox_to_anchor=(0.5, -0.03), ncol=3, frameon=True, fontsize=10)
            
            # Adjust layout
            plt.tight_layout(rect=[0, 0.08, 1, 0.96])
            
            # Save or show
            if save_path:
                save_path = Path(save_path)
                if save_path.suffix in ['.svg', '.png', '.pdf', '.jpg']:
                    save_file = save_path
                else:
                    save_path.mkdir(parents=True, exist_ok=True)
                    save_file = save_path / f"Core_{component}_all_species_COG.svg"
                
                plt.savefig(save_file, dpi=300, bbox_inches='tight')
                logger.info(f"Combined plot saved to {save_file}")
                plt.close()
            else:
                plt.show()
        
        else:
            # Original behavior: individual plots
            for species in species_with_component:
                try:
                    self.view_iModulon_weights(
                        species=species,
                        component=component,
                        save_path=save_path,
                        fig_size=fig_size,
                        font_path=font_path,
                        show_COG=show_COG
                    )
                    logger.info(f"✓ Generated plot for {species}")
                except Exception as e:
                    logger.warning(f"Failed to generate plot for {species}: {str(e)}")
        
        logger.info(f"Completed generating plots for core component '{component}'")
    
    def compare_core_iModulon(self, component: str, y_label: str = 'Species', 
                              save_path: Optional[str] = None, fig_size: Tuple[float, float] = (12, 8),
                              font_path: Optional[str] = None, reference_order: Optional[List[str]] = None) -> Dict[str, List[str]]:
        """
        Compare a core iModulon component across species with a dual-layer heatmap.
        
        Creates a heatmap showing:
        - Component membership (filled blocks from presence_matrix)
        - Gene presence across species (edge visibility from combined_gene_db)
        
        The heatmap is divided into two sections:
        - Left: genes shared by all species
        - Right: genes not shared by all species (with green edges showing presence)
        
        Parameters
        ----------
        component : str
            Core component name (e.g., 'Core_1', 'Core_2')
        y_label : str, optional
            Y-axis label. Default: 'Species'
        save_path : str, optional
            Path to save the plot. If None, displays the plot
        fig_size : tuple, optional
            Figure size as (width, height). Default: (12, 8)
        font_path : str, optional
            Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
            If provided, uses this font for all text elements
        reference_order : list of str, optional
            Custom order for species on y-axis (top to bottom).
            Example: ['MG1655', 'W', 'W3110', 'BL21', 'C', 'Crooks']
            If not provided, species are displayed in their default order.
            
        Returns
        -------
        Dict[str, List[str]]
            Dictionary with species names as keys and lists of gene names (M matrix indices) as values
            
        Raises
        ------
        ValueError
            If component is not found in any species
        """
        import matplotlib.patches as patches
        from matplotlib.collections import PatchCollection
        
        # Get all species with this component
        species_with_component = []
        result_dict = {}
        
        for species in self._species_data.keys():
            species_data = self._species_data[species]
            if species_data.M is not None and component in species_data.M.columns:
                species_with_component.append(species)
                
                # Get genes in this component for this species (using presence_matrix)
                if species_data._presence_matrix is not None:
                    component_genes = species_data._presence_matrix[
                        species_data._presence_matrix[component] == 1
                    ].index.tolist()
                    result_dict[species] = component_genes
        
        if not species_with_component:
            raise ValueError(f"Component '{component}' not found in any species")
        
        if not component.startswith('Core_'):
            warnings.warn(f"Component '{component}' does not appear to be a core component")
        
        # Apply reference_order if provided
        if reference_order:
            # Filter to only include species that have the component
            ordered_species = [sp for sp in reference_order if sp in species_with_component]
            # Add any remaining species not in reference_order
            remaining_species = [sp for sp in species_with_component if sp not in ordered_species]
            species_with_component = ordered_species + remaining_species
        
        # Set font properties if provided
        if font_path and os.path.exists(font_path):
            font_prop = fm.FontProperties(fname=font_path)
            plt.rcParams['font.family'] = font_prop.get_name()
        
        # Collect all unique genes across species for this component
        all_genes = set()
        for species in species_with_component:
            if species in result_dict:
                all_genes.update(result_dict[species])
        all_genes = sorted(list(all_genes))
        
        # Build gene presence information across species
        gene_species_presence = {}  # gene -> set of species that have this gene
        for gene in all_genes:
            gene_species_presence[gene] = set()
            
            if self.combined_gene_db is not None:
                # Check each row in combined_gene_db
                for idx, row in self.combined_gene_db.iterrows():
                    # Find leftmost gene
                    leftmost_gene = None
                    for col in self.combined_gene_db.columns:
                        val = row[col]
                        if pd.notna(val) and val != "None" and val is not None:
                            leftmost_gene = val
                            break
                    
                    # If this row contains our gene
                    if leftmost_gene == gene:
                        # Check which species have this gene
                        for species in species_with_component:
                            if species in self.combined_gene_db.columns:
                                species_gene = row[species]
                                if pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                                    gene_species_presence[gene].add(species)
                        break
        
        # Separate genes into two groups
        genes_in_all = []
        genes_not_in_all = []
        
        for gene in all_genes:
            if len(gene_species_presence.get(gene, set())) == len(species_with_component):
                genes_in_all.append(gene)
            else:
                genes_not_in_all.append(gene)
        
        # Combine with gap
        gap_width = 1
        ordered_genes = genes_in_all + genes_not_in_all
        gap_position = len(genes_in_all) if genes_in_all else -1
        
        # Create figure
        fig, ax = plt.subplots(figsize=fig_size)
        
        # Prepare data for visualization
        n_species = len(species_with_component)
        n_genes = len(ordered_genes)
        total_width = n_genes + (gap_width if gap_position > 0 else 0)
        
        # Function to create merged rectangles for adjacent species
        def create_merged_rectangles(gene_idx, x_pos, species_list, in_component_flags, exists_flags, is_right_section=False):
            rectangles = []
            i = 0
            while i < len(species_list):
                # For right section, we need to handle non-existent genes too
                if not exists_flags[i] and is_right_section:
                    # Gene doesn't exist in this species - add white rectangle for right section
                    start_i = i
                    # Find consecutive species where gene doesn't exist
                    while i < len(species_list) and not exists_flags[i]:
                        i += 1
                    
                    height = i - start_i
                    y_pos = len(species_list) - i  # Flip to have first species at top
                    
                    # White rectangle for non-existent genes
                    rect = patches.Rectangle((x_pos, y_pos), 1, height, linewidth=0,
                                           edgecolor='none', facecolor='white')
                    rectangles.append(rect)
                elif not exists_flags[i]:
                    # Left section - skip non-existent genes
                    i += 1
                    continue
                else:
                    # Gene exists - process as before
                    # Start of a group
                    start_i = i
                    start_in_component = in_component_flags[i]
                    
                    # Find consecutive species with same status
                    while i < len(species_list) and exists_flags[i] and in_component_flags[i] == start_in_component:
                        i += 1
                    
                    # Create merged rectangle
                    height = i - start_i
                    y_pos = len(species_list) - i  # Flip to have first species at top
                    
                    if start_in_component:
                        # In component - filled light blue, no edge
                        rect = patches.Rectangle((x_pos, y_pos), 1, height, linewidth=0,
                                               edgecolor='none', facecolor='#CCE5FF')  # Light blue (slightly deeper)
                    else:
                        # Not in component but exists - very light grey fill, no edge
                        rect = patches.Rectangle((x_pos, y_pos), 1, height, linewidth=0,
                                               edgecolor='none', facecolor='#F0F0F0')  # Very light grey
                    rectangles.append(rect)
            
            return rectangles
        
        # Create patches for the heatmap
        all_rectangles = []
        
        # Process each gene
        for j, gene in enumerate(ordered_genes):
            # Calculate x position accounting for gap
            if gap_position > 0 and j >= gap_position:
                x_pos = j + gap_width
            else:
                x_pos = j
            
            # Determine if this is in the right section
            is_right_section = gene in genes_not_in_all
            
            # Collect status for all species
            in_component_flags = []
            exists_flags = []
            
            for species in species_with_component:
                in_component = gene in result_dict.get(species, [])
                exists_in_species = species in gene_species_presence.get(gene, set())
                
                in_component_flags.append(in_component)
                exists_flags.append(exists_in_species)
            
            # Create merged rectangles
            rectangles = create_merged_rectangles(j, x_pos, species_with_component, 
                                                 in_component_flags, exists_flags, is_right_section)
            all_rectangles.extend(rectangles)
        
        # Add all rectangles to the plot
        for rect in all_rectangles:
            ax.add_patch(rect)
        
        # Set axis properties with small margin to show box edges
        margin = 0.1
        ax.set_xlim(-margin, total_width + margin)
        ax.set_ylim(-margin, n_species + margin)
        ax.set_aspect('equal')
        
        # No vertical line needed - gap will separate the sections
        
        # Set ticks and labels
        x_positions = []
        for j in range(len(ordered_genes)):
            if gap_position > 0 and j >= gap_position:
                x_positions.append(j + gap_width + 0.5)
            else:
                x_positions.append(j + 0.5)
        
        ax.set_xticks(x_positions)
        ax.set_yticks(np.arange(n_species) + 0.5)
        
        # Get gene names for x-axis labels
        x_labels = []
        for gene in ordered_genes:
            # Try to get gene_name from gene_table
            gene_name = gene  # Default to using the index
            for species in species_with_component:
                species_data = self._species_data[species]
                if species_data.gene_table is not None:
                    # Map from M matrix index to species gene
                    if self.combined_gene_db is not None and species in self.combined_gene_db.columns:
                        # Find the species-specific gene name
                        for idx, row in self.combined_gene_db.iterrows():
                            leftmost_gene = None
                            for col in self.combined_gene_db.columns:
                                val = row[col]
                                if pd.notna(val) and val != "None" and val is not None:
                                    leftmost_gene = val
                                    break
                            
                            if leftmost_gene == gene:
                                species_gene = row[species]
                                if pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                                    if species_gene in species_data.gene_table.index:
                                        if 'gene_name' in species_data.gene_table.columns:
                                            name = species_data.gene_table.loc[species_gene, 'gene_name']
                                            if pd.notna(name) and name != '':
                                                gene_name = name
                                                break
                    
                    # Also check if gene exists directly in gene_table
                    elif gene in species_data.gene_table.index:
                        if 'gene_name' in species_data.gene_table.columns:
                            name = species_data.gene_table.loc[gene, 'gene_name']
                            if pd.notna(name) and name != '':
                                gene_name = name
                                break
            
            x_labels.append(gene_name)
        
        ax.set_xticklabels(x_labels, rotation=45, ha='right')
        # Reverse y-axis labels to match top-to-bottom order
        ax.set_yticklabels(species_with_component[::-1])
        
        # Set labels and title
        ax.set_xlabel('Genes', fontsize=12)
        ax.set_ylabel(y_label, fontsize=12)
        ax.set_title(f'Core iModulon {component} Comparison', fontsize=14)
        
        # Apply font to all text elements if font_path provided
        if font_path and os.path.exists(font_path):
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontproperties(font_prop)
            ax.xaxis.label.set_fontproperties(font_prop)
            ax.yaxis.label.set_fontproperties(font_prop)
            ax.title.set_fontproperties(font_prop)
        
        # Remove all spines first
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        # Add separate boxes for left and right heatmaps
        if gap_position > 0:
            # Left heatmap box
            left_box = patches.Rectangle((0, 0), gap_position, n_species, 
                                       linewidth=1, edgecolor='black', facecolor='none')
            ax.add_patch(left_box)
            
            # Right heatmap box
            right_start = gap_position + gap_width
            right_width = total_width - right_start
            right_box = patches.Rectangle((right_start, 0), right_width, n_species,
                                        linewidth=1, edgecolor='black', facecolor='none')
            ax.add_patch(right_box)
        else:
            # Single box if no separation
            full_box = patches.Rectangle((0, 0), total_width, n_species,
                                       linewidth=1, edgecolor='black', facecolor='none')
            ax.add_patch(full_box)
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#CCE5FF', edgecolor='black', linewidth=1, label='Present in iModulon'),
            Patch(facecolor='#F0F0F0', edgecolor='black', linewidth=1, label='Absent from iModulon'),
            Patch(facecolor='white', edgecolor='black', linewidth=1, label=f'Absent from {y_label.lower()}')
        ]
        
        # Position legend to the right with same gap as between heatmaps
        if gap_position > 0:
            legend_x = total_width + gap_width
        else:
            legend_x = total_width + 1
            
        # Create custom legend
        legend = ax.legend(handles=legend_elements, 
                          loc='center left',
                          bbox_to_anchor=(legend_x / total_width, 0.5),
                          frameon=True,
                          framealpha=1,
                          edgecolor='black')
        
        # Apply font to legend if provided
        if font_path and os.path.exists(font_path):
            for text in legend.get_texts():
                text.set_fontproperties(font_prop)
        
        # Adjust figure size to accommodate legend
        fig.subplots_adjust(right=0.75)
        
        # Tight layout
        plt.tight_layout()
        
        # Save or show
        if save_path:
            save_path = Path(save_path)
            if save_path.suffix in ['.svg', '.png', '.pdf', '.jpg']:
                save_file = save_path
            else:
                save_path.mkdir(parents=True, exist_ok=True)
                save_file = save_path / f"{component}_comparison_heatmap.svg"
            
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            logger.info(f"Comparison heatmap saved to {save_file}")
            
        # Always show the plot (useful for Jupyter notebooks)
        plt.show()
        
        return result_dict
    
    def save_to_json_multimodulon(self, save_path: str):
        """
        Save the entire MultiModulon object to a JSON file.
        
        This method serializes all data including:
        - All species data (log_tpm, log_tpm_norm, X, M, A, gene_table, sample_sheet, M_thresholds, presence_matrix)
        - Combined gene database
        - Input folder path
        
        Parameters
        ----------
        save_path : str
            Path to save the JSON file. If the path ends with '.json.gz', 
            the file will be gzip compressed.
            
        Examples
        --------
        >>> multiModulon.save_to_json_multimodulon("multimodulon_data.json")
        >>> multiModulon.save_to_json_multimodulon("multimodulon_data.json.gz")  # Compressed
        """
        import gzip
        
        save_path = Path(save_path)
        
        # Create the data structure to save
        data = {
            'input_folder_path': str(self.input_folder_path),
            'species_data': {},
            'combined_gene_db': None
        }
        
        # Save combined gene database if available
        if hasattr(self, '_combined_gene_db') and self._combined_gene_db is not None:
            data['combined_gene_db'] = {
                'data': self._combined_gene_db.to_dict('records'),
                'columns': list(self._combined_gene_db.columns)
            }
        
        
        # Save each species data
        for species_name, species_data_obj in self._species_data.items():
            species_dict = {
                'species_name': species_data_obj.species_name,
                'data_path': str(species_data_obj.data_path),
                'log_tpm': None,
                'log_tpm_norm': None,
                'X': None,
                'M': None,
                'A': None,
                'sample_sheet': None,
                'gene_table': None,
                'M_thresholds': None,
                'presence_matrix': None
            }
            
            # Save DataFrames if they exist
            if species_data_obj._log_tpm is not None:
                species_dict['log_tpm'] = {
                    'data': species_data_obj._log_tpm.to_dict('records'),
                    'index': list(species_data_obj._log_tpm.index),
                    'columns': list(species_data_obj._log_tpm.columns)
                }
            
            if species_data_obj._log_tpm_norm is not None:
                species_dict['log_tpm_norm'] = {
                    'data': species_data_obj._log_tpm_norm.to_dict('records'),
                    'index': list(species_data_obj._log_tpm_norm.index),
                    'columns': list(species_data_obj._log_tpm_norm.columns)
                }
            
            if species_data_obj._X is not None:
                species_dict['X'] = {
                    'data': species_data_obj._X.to_dict('records'),
                    'index': list(species_data_obj._X.index),
                    'columns': list(species_data_obj._X.columns)
                }
            
            if species_data_obj._M is not None:
                species_dict['M'] = {
                    'data': species_data_obj._M.to_dict('records'),
                    'index': list(species_data_obj._M.index),
                    'columns': list(species_data_obj._M.columns)
                }
            
            if species_data_obj._A is not None:
                species_dict['A'] = {
                    'data': species_data_obj._A.to_dict('records'),
                    'index': list(species_data_obj._A.index),
                    'columns': list(species_data_obj._A.columns)
                }
            
            if species_data_obj._sample_sheet is not None:
                species_dict['sample_sheet'] = {
                    'data': species_data_obj._sample_sheet.to_dict('records'),
                    'index': list(species_data_obj._sample_sheet.index),
                    'columns': list(species_data_obj._sample_sheet.columns)
                }
            
            if species_data_obj._gene_table is not None:
                species_dict['gene_table'] = {
                    'data': species_data_obj._gene_table.to_dict('records'),
                    'index': list(species_data_obj._gene_table.index),
                    'columns': list(species_data_obj._gene_table.columns)
                }
            
            if species_data_obj._M_thresholds is not None:
                species_dict['M_thresholds'] = {
                    'data': species_data_obj._M_thresholds.to_dict('records'),
                    'index': list(species_data_obj._M_thresholds.index),
                    'columns': list(species_data_obj._M_thresholds.columns)
                }
            
            if species_data_obj._presence_matrix is not None:
                species_dict['presence_matrix'] = {
                    'data': species_data_obj._presence_matrix.to_dict('records'),
                    'index': list(species_data_obj._presence_matrix.index),
                    'columns': list(species_data_obj._presence_matrix.columns)
                }
            
            data['species_data'][species_name] = species_dict
        
        # Save to JSON file with optional compression
        if str(save_path).endswith('.json.gz'):
            # Save as compressed gzip file
            with gzip.open(save_path, 'wt', encoding='utf-8') as f:
                json.dump(data, f, indent=2)
            logger.info(f"MultiModulon object saved to {save_path} (compressed)")
        else:
            # Save as regular JSON file
            with open(save_path, 'w') as f:
                json.dump(data, f, indent=2)
            logger.info(f"MultiModulon object saved to {save_path}")
    
    @staticmethod
    def load_json_multimodulon(load_path: str) -> 'MultiModulon':
        """
        Load a MultiModulon object from a JSON file.
        
        This static method creates a new MultiModulon object and populates it
        with all the data saved in the JSON file.
        
        Parameters
        ----------
        load_path : str
            Path to the JSON file to load. If the path ends with '.json.gz',
            the file will be treated as gzip compressed.
            
        Returns
        -------
        MultiModulon
            A new MultiModulon object with all loaded data
            
        Examples
        --------
        >>> multiModulon = MultiModulon.load_json_multimodulon("multimodulon_data.json")
        >>> multiModulon = MultiModulon.load_json_multimodulon("multimodulon_data.json.gz")  # Compressed
        """
        import gzip
        
        load_path = Path(load_path)
        
        # Load JSON data with optional decompression
        if str(load_path).endswith('.json.gz'):
            # Load from compressed gzip file
            with gzip.open(load_path, 'rt', encoding='utf-8') as f:
                data = json.load(f)
        else:
            # Load from regular JSON file
            with open(load_path, 'r') as f:
                data = json.load(f)
        
        # Create a new MultiModulon object
        # We need to handle the case where the input folder might not exist
        input_folder = data['input_folder_path']
        
        # Create a minimal MultiModulon object
        multi_modulon = MultiModulon.__new__(MultiModulon)
        multi_modulon.input_folder_path = Path(input_folder)
        multi_modulon._species_data = {}
        multi_modulon._bbh = None
        
        # Load combined gene database
        if data.get('combined_gene_db') is not None:
            df_data = data['combined_gene_db']
            multi_modulon._combined_gene_db = pd.DataFrame(df_data['data'])
            if df_data['data']:  # Only set columns if there's data
                multi_modulon._combined_gene_db = multi_modulon._combined_gene_db[df_data['columns']]
        else:
            multi_modulon._combined_gene_db = None
        
        
        # Load species data
        for species_name, species_dict in data['species_data'].items():
            # Create SpeciesData object
            species_data = SpeciesData.__new__(SpeciesData)
            species_data.species_name = species_dict['species_name']
            species_data.data_path = Path(species_dict['data_path'])
            
            # Load DataFrames
            def load_dataframe(df_data):
                if df_data is None:
                    return None
                df = pd.DataFrame(df_data['data'])
                if df_data['data']:  # Only set index/columns if there's data
                    df.index = df_data['index']
                    df.columns = df_data['columns']
                return df
            
            species_data._log_tpm = load_dataframe(species_dict.get('log_tpm'))
            species_data._log_tpm_norm = load_dataframe(species_dict.get('log_tpm_norm'))
            species_data._X = load_dataframe(species_dict.get('X'))
            species_data._M = load_dataframe(species_dict.get('M'))
            species_data._A = load_dataframe(species_dict.get('A'))
            species_data._sample_sheet = load_dataframe(species_dict.get('sample_sheet'))
            species_data._gene_table = load_dataframe(species_dict.get('gene_table'))
            species_data._M_thresholds = load_dataframe(species_dict.get('M_thresholds'))
            species_data._presence_matrix = load_dataframe(species_dict.get('presence_matrix'))
            
            # Add to MultiModulon
            multi_modulon._species_data[species_name] = species_data
        
        logger.info(f"MultiModulon object loaded from {load_path}")
        logger.info(f"Loaded {len(multi_modulon._species_data)} species")
        
        return multi_modulon