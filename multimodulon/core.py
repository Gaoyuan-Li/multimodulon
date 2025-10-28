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
from .optimization import optimize_number_of_core_components, optimize_number_of_unique_components, passes_single_gene_filter
from .utils import BBHAnalyzer
from .gff_utils import gff2pandas, create_gene_table, extract_protein_sequences
from .multiview_ica import run_multiview_ica
from .core_io import save_bbh, load_bbh, get_orthologs, save_to_json_multimodulon, load_json_multimodulon
from .plotting import (
    view_iModulon_weights,
    view_iModulon_activities,
    view_iModulon_genes,
    view_core_iModulon_weights,
    compare_core_iModulon,
    compare_core_iModulon_activity,
    show_iModulon_activity_change,
    show_gene_iModulon_correlation,
    plot_iM_conservation_bubble_matrix,
)
from sklearn.cluster import HDBSCAN, DBSCAN
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch
import os
from scipy.stats import pearsonr, median_abs_deviation
from scipy.signal import argrelextrema

logger = logging.getLogger(__name__)


def generate_species_colors(species_list):
    """
    Generate visually distinct colors for a list of species.
    
    Uses HSL color space to create colors with good visual separation,
    low saturation (30%), and moderate lightness (50%).
    
    Parameters
    ----------
    species_list : list
        List of species names
        
    Returns
    -------
    dict
        Dictionary mapping species names to hex color codes
    """
    import colorsys
    
    n_species = len(species_list)
    colors = {}
    
    # Use golden ratio for hue distribution to maximize visual separation
    golden_ratio = 0.618033988749895
    hue = 0
    
    for i, species in enumerate(sorted(species_list)):
        # Convert HSL to RGB (hue in [0,1], saturation=0.4, lightness=0.5)
        rgb = colorsys.hls_to_rgb(hue, 0.5, 0.4)
        # Convert to hex
        hex_color = '#{:02X}{:02X}{:02X}'.format(
            int(rgb[0] * 255),
            int(rgb[1] * 255),
            int(rgb[2] * 255)
        )
        colors[species] = hex_color
        
        # Increment hue using golden ratio
        hue = (hue + golden_ratio) % 1.0
    
    return colors


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
        
        # Initialize species color palette
        self.species_palette = generate_species_colors(list(self._species_data.keys()))
    
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
        
        # Prepare X matrices dictionary
        species_X_matrices = {}
        for species in species_list:
            species_X_matrices[species] = self._species_data[species].X
        
        # Run multi-view ICA
        results = run_multiview_ica(
            species_X_matrices=species_X_matrices,
            a_values=a_values,
            c=c,
            mode=mode
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
        seed: int = 42
    ) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
        """
        Run robust multi-view ICA with clustering to identify consistent components.
        
        This method runs multi-view ICA multiple times and uses clustering to identify
        robust components that appear consistently across runs. Components whose largest
        absolute weight is less than three times the second-largest (per species) are
        treated as single-gene components and removed before clustering.
        
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
        Components are filtered using the single-gene criterion described in the
        package overview; no additional thresholds are required.
            
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
        ...     num_runs=100
        ... )
        
        >>> # Access the results
        >>> M_strain1 = M_matrices['strain1']  # Mixing matrix for strain1
        >>> A_strain1 = A_matrices['strain1']  # Activity matrix for strain1
        """
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
                    
                    if passes_single_gene_filter(weight_vector):
                        component = self._enforce_sign_convention(weight_vector)
                        core_components[species].append(component)
                
                # Process unique components  
                for comp_idx in range(c, a[species]):
                    weight_vector = M.iloc[:, comp_idx].values
                    
                    if passes_single_gene_filter(weight_vector):
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
                    component_vector = np.asarray(centroids[species])
                    if passes_single_gene_filter(component_vector):
                        components.append(component_vector)
                        component_names.append(core_name)
            
            # Add unique components
            for unique_name, centroid in unique_results[species]:
                component_vector = np.asarray(centroid)
                if passes_single_gene_filter(component_vector):
                    components.append(component_vector)
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
        any_species = species_list[0]
        core_count = len([col for col in final_M_matrices[any_species].columns if col.startswith('Core_')])
        print(f"Total core components retained: {core_count}")
        for species in species_list:
            n_unique = len([col for col in final_M_matrices[species].columns if col.startswith('Unique_')])
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
    
    def optimize_number_of_core_components(self, **kwargs) -> int:
        """Optimize the number of core components using the single-gene filter."""
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
    
    def save_gene_table(self, output_folder: str) -> None:
        """
        Save gene tables for all species to CSV files.
        
        This method saves the gene table of each species to a separate CSV file
        in the specified folder. The files are named using the species name.
        
        Parameters
        ----------
        output_folder : str
            Path to the folder where gene tables will be saved.
            The folder will be created if it doesn't exist.
            
        Notes
        -----
        Each gene table is saved as a CSV file with the following naming pattern:
        - {species_name}.csv
        
        The CSV files include the locus_tag as the index column.
        """
        output_path = Path(output_folder)
        
        # Create output folder if it doesn't exist
        if not output_path.exists():
            logger.info(f"Creating output folder: {output_path}")
            output_path.mkdir(parents=True, exist_ok=True)
            print(f"Created output folder: {output_path}")
        
        print(f"\nSaving gene tables to {output_folder}")
        print("=" * 60)
        
        saved_count = 0
        
        for species_name, species_data in self._species_data.items():
            print(f"\nProcessing {species_name}...")
            
            # Check if gene_table exists
            if species_data._gene_table is None:
                logger.warning(f"Gene table not found for {species_name}. Skipping...")
                print(f"  ✗ Gene table not found. Skipping...")
                continue
            
            # Define output file path
            output_file = output_path / f"{species_name}.csv"
            
            try:
                # Save gene table to CSV
                species_data._gene_table.to_csv(output_file)
                saved_count += 1
                
                # Get table info
                num_genes = len(species_data._gene_table)
                num_columns = len(species_data._gene_table.columns)
                
                print(f"  ✓ Saved {num_genes} genes × {num_columns} columns to {output_file.name}")
                logger.info(f"Saved gene table for {species_name}: {num_genes} genes, {num_columns} columns")
                
            except Exception as e:
                logger.error(f"Error saving gene table for {species_name}: {e}")
                print(f"  ✗ Error: {e}")
        
        print("\n" + "=" * 60)
        print(f"Gene table saving completed! Saved {saved_count} species tables.")
    
    
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
            
            # Create presence matrix with mapped indexes
            # First, create a mapping from M matrix indexes (leftmost genes) to species-specific genes
            gene_mapping = {}
            
            if self.combined_gene_db is not None and species_name in self.combined_gene_db.columns:
                for _, row in self.combined_gene_db.iterrows():
                    # Find leftmost gene
                    leftmost_gene = None
                    for col in self.combined_gene_db.columns:
                        val = row[col]
                        if pd.notna(val) and val != "None" and val is not None:
                            leftmost_gene = val
                            break
                    
                    # Get species-specific gene
                    if leftmost_gene and leftmost_gene in M.index:
                        species_gene = row[species_name]
                        if pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                            # Only include if gene exists in gene_table
                            if species_data.gene_table is not None and species_gene in species_data.gene_table.index:
                                gene_mapping[leftmost_gene] = species_gene
            else:
                # No combined_gene_db, use genes directly if they exist in gene_table
                if species_data.gene_table is not None:
                    for gene in M.index:
                        if gene in species_data.gene_table.index:
                            gene_mapping[gene] = gene
            
            # Create presence matrix with ALL genes from gene_table
            if species_data.gene_table is not None:
                # Initialize presence matrix with all genes from gene_table, all values as 0
                presence_matrix = pd.DataFrame(
                    0,  # Initialize all values to 0
                    index=species_data.gene_table.index,
                    columns=M.columns,
                    dtype=int
                )
                
                # Binarize based on thresholds using the mapping
                for component in M.columns:
                    threshold = thresholds[component]
                    # Map and binarize only for genes that have mapping
                    for leftmost_gene, species_gene in gene_mapping.items():
                        if leftmost_gene in M.index and species_gene in presence_matrix.index:
                            presence_matrix.loc[species_gene, component] = int(abs(M.loc[leftmost_gene, component]) > threshold)
                
                # Report statistics
                n_mapped_genes = len(gene_mapping)
                n_total_genes = len(presence_matrix)
                print(f"  Gene mapping: {n_mapped_genes}/{n_total_genes} genes have expression data")
            else:
                # Fallback to original behavior if no gene_table
                presence_matrix = pd.DataFrame(
                    index=list(gene_mapping.values()),
                    columns=M.columns,
                    dtype=int
                )
                
                for component in M.columns:
                    threshold = thresholds[component]
                    for leftmost_gene, species_gene in gene_mapping.items():
                        if leftmost_gene in M.index:
                            presence_matrix.loc[species_gene, component] = int(abs(M.loc[leftmost_gene, component]) > threshold)
            
            # Store results
            species_data._M_thresholds = M_thresholds_df
            species_data._presence_matrix = presence_matrix
            
            # Report statistics
            n_components = len(M.columns)
            avg_genes_per_component = presence_matrix.sum().mean()
            print(f"✓ Optimized thresholds for {n_components} components")
            print(f"  Average genes per component: {avg_genes_per_component:.1f}")
        
        print("\nThreshold optimization completed!")
    
    def update_M_threshold(self, component: str, species: str, new_threshold: float):
        """
        Update the M threshold for a specific component in a specific species.
        
        This method updates both the M_thresholds and presence_matrix for the
        specified component using the new threshold value.
        
        Parameters
        ----------
        component : str
            Name of the component to update
        species : str
            Name of the species to update
        new_threshold : float
            New threshold value to apply
            
        Raises
        ------
        ValueError
            If species not found, component not found, or M/M_thresholds not available
        """
        # Validate species
        if species not in self._species_data:
            raise ValueError(f"Species '{species}' not found in MultiModulon object")
        
        species_data = self._species_data[species]
        
        # Check if M matrix is available
        if not hasattr(species_data, '_M') or species_data._M is None:
            raise ValueError(f"M matrix not available for species '{species}'. Run compute_M() first.")
        
        # Check if M_thresholds is available
        if not hasattr(species_data, '_M_thresholds') or species_data._M_thresholds is None:
            raise ValueError(f"M_thresholds not available for species '{species}'. Run optimize_M_thresholds() first.")
        
        # Check if component exists
        M = species_data.M
        if component not in M.columns:
            raise ValueError(f"Component '{component}' not found in M matrix for species '{species}'")
        
        # Update M_thresholds with explicit dtype casting to match DataFrame dtype
        species_data._M_thresholds.loc[component, 'M_threshold'] = np.float32(new_threshold)
        
        # Update presence_matrix
        if hasattr(species_data, '_presence_matrix') and species_data._presence_matrix is not None:
            # Need to use the same gene mapping as in optimize_M_thresholds
            # Create mapping from M matrix indexes (leftmost genes) to species-specific genes
            gene_mapping = {}
            
            if self.combined_gene_db is not None and species in self.combined_gene_db.columns:
                for _, row in self.combined_gene_db.iterrows():
                    # Find leftmost gene
                    leftmost_gene = None
                    for col in self.combined_gene_db.columns:
                        val = row[col]
                        if pd.notna(val) and val != "None" and val is not None:
                            leftmost_gene = val
                            break
                    
                    # Get species-specific gene
                    if leftmost_gene and leftmost_gene in M.index:
                        species_gene = row[species]
                        if pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                            # Only include if gene exists in gene_table
                            if species_data.gene_table is not None and species_gene in species_data.gene_table.index:
                                gene_mapping[leftmost_gene] = species_gene
            else:
                # No combined_gene_db, use genes directly if they exist in gene_table
                if species_data.gene_table is not None:
                    for gene in M.index:
                        if gene in species_data.gene_table.index:
                            gene_mapping[gene] = gene
            
            # First reset all values to 0 for this component
            species_data._presence_matrix[component] = 0
            
            # Then set to 1 for genes above threshold using the mapping
            for leftmost_gene, species_gene in gene_mapping.items():
                if species_gene in species_data._presence_matrix.index and leftmost_gene in M.index:
                    species_data._presence_matrix.loc[species_gene, component] = int(abs(M.loc[leftmost_gene, component]) > new_threshold)
            
            # Report the change
            n_genes = species_data._presence_matrix[component].sum()
            print(f"Updated threshold for component '{component}' in species '{species}':")
            print(f"  New threshold: {new_threshold:.4f}")
            print(f"  Number of genes above threshold: {n_genes}")
        else:
            print(f"Warning: presence_matrix not found for species '{species}'. Only M_thresholds updated.")
    
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
    
    
    @staticmethod

    # Wrapper methods for IO operations
    def save_bbh(self, output_path: str):
        """Save BBH results to file."""
        return save_bbh(self, output_path)
    
    def load_bbh(self, input_path: str):
        """Load BBH results from file."""
        return load_bbh(self, input_path)
    
    def get_orthologs(self, species1: str, species2: str):
        """Get ortholog pairs between two species."""
        return get_orthologs(self, species1, species2)
    
    def rename_core_iModulon(self, old_name: str, new_name: str):
        """
        Rename a core iModulon across all species and related dataframes.
        
        Parameters
        ----------
        old_name : str
            Current name of the core iModulon (e.g., 'Core_1')
        new_name : str
            New name for the core iModulon (e.g., 'ArgR')
            
        Raises
        ------
        ValueError
            If old_name doesn't exist in any species or new_name already exists
        """
        # Check if the core iModulon exists in any species
        found_in_species = []
        for species_name in self.species:
            species_data = self[species_name]
            if species_data.M is not None and old_name in species_data.M.columns:
                found_in_species.append(species_name)
        
        if not found_in_species:
            raise ValueError(f"Core iModulon '{old_name}' not found in any species")
        
        # Check if new name already exists in any species
        for species_name in self.species:
            species_data = self[species_name]
            if species_data.M is not None and new_name in species_data.M.columns:
                raise ValueError(f"iModulon '{new_name}' already exists in species '{species_name}'")
        
        # Rename in all species
        for species_name in self.species:
            species_data = self[species_name]
            
            # Skip if this species doesn't have the iModulon
            if species_data.M is None or old_name not in species_data.M.columns:
                print(f"Warning: Core iModulon '{old_name}' not found in species '{species_name}'")
                continue
            
            # Rename in M matrix (columns)
            if species_data._M is not None:
                species_data._M = species_data._M.rename(columns={old_name: new_name})
            
            # Rename in A matrix (rows)
            if species_data._A is not None:
                species_data._A = species_data._A.rename(index={old_name: new_name})
            
            # Rename in M_thresholds (index)
            if hasattr(species_data, '_M_thresholds') and species_data._M_thresholds is not None:
                species_data._M_thresholds = species_data._M_thresholds.rename(index={old_name: new_name})
            
            # Rename in presence_matrix (columns)
            if hasattr(species_data, '_presence_matrix') and species_data._presence_matrix is not None:
                species_data._presence_matrix = species_data._presence_matrix.rename(columns={old_name: new_name})
        
        print(f"Successfully renamed core iModulon '{old_name}' to '{new_name}' in {len(found_in_species)} species")
    
    def rename_unique_iModulon(self, species: str, old_name: str, new_name: str):
        """
        Rename a unique iModulon for a specific species.
        
        Parameters
        ----------
        species : str
            Species name where the unique iModulon exists
        old_name : str
            Current name of the unique iModulon (e.g., 'Unique_1')
        new_name : str
            New name for the unique iModulon
            
        Raises
        ------
        ValueError
            If species doesn't exist, old_name doesn't exist in the species, or new_name already exists
        """
        # Validate species
        if species not in self.species:
            raise ValueError(f"Species '{species}' not found. Available species: {', '.join(self.species)}")
        
        # Get species data
        species_data = self[species]
        
        # Check if the unique iModulon exists
        if species_data.M is None or old_name not in species_data.M.columns:
            raise ValueError(f"iModulon '{old_name}' not found in species '{species}'")
        
        # Check if new name already exists
        if species_data.M is not None and new_name in species_data.M.columns:
            raise ValueError(f"iModulon '{new_name}' already exists in species '{species}'")
        
        # Rename in M matrix (columns)
        if species_data._M is not None:
            species_data._M = species_data._M.rename(columns={old_name: new_name})
        
        # Rename in A matrix (rows)
        if species_data._A is not None:
            species_data._A = species_data._A.rename(index={old_name: new_name})
        
        # Rename in M_thresholds (index)
        if hasattr(species_data, '_M_thresholds') and species_data._M_thresholds is not None:
            species_data._M_thresholds = species_data._M_thresholds.rename(index={old_name: new_name})
        
        # Rename in presence_matrix (columns)
        if hasattr(species_data, '_presence_matrix') and species_data._presence_matrix is not None:
            species_data._presence_matrix = species_data._presence_matrix.rename(columns={old_name: new_name})
        
        print(f"Successfully renamed unique iModulon '{old_name}' to '{new_name}' in species '{species}'")
    
    def save_to_json_multimodulon(self, save_path: str):
        """Save MultiModulon object to JSON format."""
        return save_to_json_multimodulon(self, save_path)
    
    @staticmethod
    def load_json_multimodulon(load_path: str):
        """Load MultiModulon object from JSON format."""
        return load_json_multimodulon(load_path)
    
    # Wrapper methods for plotting operations
    def view_iModulon_weights(self, *args, **kwargs):
        """View weights of an iModulon."""
        return view_iModulon_weights(self, *args, **kwargs)
    
    def view_iModulon_activities(self, *args, **kwargs):
        """View activities of an iModulon."""
        return view_iModulon_activities(self, *args, **kwargs)
    
    def view_iModulon_genes(self, *args, **kwargs):
        """View genes in an iModulon."""
        return view_iModulon_genes(self, *args, **kwargs)
    
    def view_core_iModulon_weights(self, *args, **kwargs):
        """View weights of a core iModulon."""
        return view_core_iModulon_weights(self, *args, **kwargs)
    
    def compare_core_iModulon(self, *args, **kwargs):
        """Compare core iModulon across species."""
        return compare_core_iModulon(self, *args, **kwargs)
    
    def plot_iM_conservation_bubble_matrix(self, *args, **kwargs):
        """Plot the cross-species iModulon conservation bubble matrix."""
        return plot_iM_conservation_bubble_matrix(self, *args, **kwargs)
    
    def compare_core_iModulon_activity(self, *args, **kwargs):
        """Compare core iModulon activities across species."""
        return compare_core_iModulon_activity(self, *args, **kwargs)
    
    def show_iModulon_activity_change(self, *args, **kwargs):
        """Visualize iModulon activity changes between two conditions."""
        return show_iModulon_activity_change(self, *args, **kwargs)
    
    def show_gene_iModulon_correlation(self, *args, **kwargs):
        """Show correlation between gene expression and iModulon activity across species."""
        return show_gene_iModulon_correlation(self, *args, **kwargs)
    
    def core_iModulon_stability(self, *args, **kwargs):
        """Quantify core iModulon stability across species."""
        from .stability import core_iModulon_stability
        return core_iModulon_stability(self, *args, **kwargs)
    
    def remove_component(self, component: str, species: Optional[str] = None):
        """
        Remove a component from the MultiModulon object and rename subsequent components
        to maintain continuous numbering.
        
        For Core components (starting with 'Core_'), removes from all species by default
        and renames Core_8 to Core_7, Core_9 to Core_8, etc.
        
        For Unique components (starting with 'Unique_'), requires species specification
        and renames subsequent Unique components in that species only.
        
        Parameters
        ----------
        component : str
            Name of the component to remove (e.g., 'Core_1', 'Unique_1')
        species : str, optional
            Species name. Required for Unique components, ignored with warning for Core components.
            
        Raises
        ------
        ValueError
            If component not found, or if species required but not provided for Unique components.
        """
        if not component:
            raise ValueError("Component name cannot be empty")
        
        # Determine component type
        is_core = component.startswith('Core_')
        is_unique = component.startswith('Unique_')
        
        if not is_core and not is_unique:
            warnings.warn(f"Component '{component}' does not follow standard naming (Core_* or Unique_*)")
        
        # Handle Core components
        if is_core:
            if species is not None:
                warnings.warn(f"Species parameter ignored for Core component '{component}'. "
                            f"Core components will be removed from all species.")
            
            # Remove from all species
            removed_from = []
            for species_name, species_data in self._species_data.items():
                if self._remove_component_from_species(component, species_name):
                    removed_from.append(species_name)
            
            if not removed_from:
                raise ValueError(f"Component '{component}' not found in any species")
            
            logger.info(f"Removed Core component '{component}' from species: {', '.join(removed_from)}")
            
            # Rename subsequent Core components in all species
            component_num = int(component.split('_')[1])
            self._rename_subsequent_components('Core', component_num, species_list=removed_from)
        
        # Handle Unique components
        elif is_unique:
            if species is None:
                raise ValueError(f"Species parameter is required for Unique component '{component}'")
            
            if species not in self._species_data:
                raise ValueError(f"Species '{species}' not found in loaded data")
            
            if not self._remove_component_from_species(component, species):
                raise ValueError(f"Component '{component}' not found in species '{species}'")
            
            logger.info(f"Removed Unique component '{component}' from species '{species}'")
            
            # Rename subsequent Unique components in this species only
            component_num = int(component.split('_')[1])
            self._rename_subsequent_components('Unique', component_num, species_list=[species])
        
        # Handle other components
        else:
            if species is None:
                # Try to remove from all species
                removed_from = []
                for species_name in self._species_data:
                    if self._remove_component_from_species(component, species_name):
                        removed_from.append(species_name)
                
                if not removed_from:
                    raise ValueError(f"Component '{component}' not found in any species")
                
                logger.info(f"Removed component '{component}' from species: {', '.join(removed_from)}")
            else:
                # Remove from specific species
                if species not in self._species_data:
                    raise ValueError(f"Species '{species}' not found in loaded data")
                
                if not self._remove_component_from_species(component, species):
                    raise ValueError(f"Component '{component}' not found in species '{species}'")
                
                logger.info(f"Removed component '{component}' from species '{species}'")
    
    def _remove_component_from_species(self, component: str, species: str) -> bool:
        """
        Helper method to remove component from a specific species.
        
        Parameters
        ----------
        component : str
            Component name to remove
        species : str
            Species name
            
        Returns
        -------
        bool
            True if component was found and removed, False otherwise
        """
        species_data = self._species_data[species]
        component_found = False
        
        # Remove from M matrix (column)
        if species_data._M is not None and component in species_data._M.columns:
            species_data._M = species_data._M.drop(columns=[component])
            component_found = True
            logger.debug(f"Removed '{component}' from M matrix of {species}")
        
        # Remove from A matrix (row)
        if species_data._A is not None and component in species_data._A.index:
            species_data._A = species_data._A.drop(index=[component])
            component_found = True
            logger.debug(f"Removed '{component}' from A matrix of {species}")
        
        # Remove from M_thresholds (row)
        if species_data._M_thresholds is not None and component in species_data._M_thresholds.index:
            species_data._M_thresholds = species_data._M_thresholds.drop(index=[component])
            component_found = True
            logger.debug(f"Removed '{component}' from M_thresholds of {species}")
        
        # Remove from presence_matrix (column)
        if species_data._presence_matrix is not None and component in species_data._presence_matrix.columns:
            species_data._presence_matrix = species_data._presence_matrix.drop(columns=[component])
            component_found = True
            logger.debug(f"Removed '{component}' from presence_matrix of {species}")
        
        return component_found
    
    def _rename_subsequent_components(self, component_type: str, removed_num: int, species_list: List[str]):
        """
        Rename subsequent components after removal to maintain continuous numbering.
        
        Parameters
        ----------
        component_type : str
            'Core' or 'Unique'
        removed_num : int
            The number of the removed component
        species_list : List[str]
            List of species to rename components in
        """
        logger.info(f"Renaming subsequent {component_type} components after removing {component_type}_{removed_num}")
        
        # Find all components that need renaming
        for species_name in species_list:
            species_data = self._species_data[species_name]
            components_to_rename = []
            
            # Check M matrix for components that need renaming
            if species_data._M is not None:
                for col in species_data._M.columns:
                    # Only process components that match the exact pattern (Core_N or Unique_N)
                    if col.startswith(f'{component_type}_'):
                        parts = col.split('_')
                        # Ensure it's exactly in format "Core_N" or "Unique_N" (not renamed components)
                        if len(parts) == 2 and parts[0] == component_type:
                            try:
                                comp_num = int(parts[1])
                                if comp_num > removed_num:
                                    components_to_rename.append((col, comp_num))
                            except ValueError:
                                # Skip if the part after underscore is not a number
                                continue
            
            # Sort components by number
            components_to_rename.sort(key=lambda x: x[1])
            
            # Rename each component
            for old_name, old_num in components_to_rename:
                new_num = old_num - 1
                new_name = f'{component_type}_{new_num}'
                
                # Rename in M matrix (column)
                if species_data._M is not None and old_name in species_data._M.columns:
                    species_data._M = species_data._M.rename(columns={old_name: new_name})
                    logger.debug(f"Renamed '{old_name}' to '{new_name}' in M matrix of {species_name}")
                
                # Rename in A matrix (row)
                if species_data._A is not None and old_name in species_data._A.index:
                    species_data._A = species_data._A.rename(index={old_name: new_name})
                    logger.debug(f"Renamed '{old_name}' to '{new_name}' in A matrix of {species_name}")
                
                # Rename in M_thresholds (row)
                if species_data._M_thresholds is not None and old_name in species_data._M_thresholds.index:
                    species_data._M_thresholds = species_data._M_thresholds.rename(index={old_name: new_name})
                    logger.debug(f"Renamed '{old_name}' to '{new_name}' in M_thresholds of {species_name}")
                
                # Rename in presence_matrix (column)
                if species_data._presence_matrix is not None and old_name in species_data._presence_matrix.columns:
                    species_data._presence_matrix = species_data._presence_matrix.rename(columns={old_name: new_name})
                    logger.debug(f"Renamed '{old_name}' to '{new_name}' in presence_matrix of {species_name}")
        
        logger.info(f"Completed renaming of {component_type} components")
