"""Core MultiModulon class for multi-species expression analysis."""

from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import json
import pickle
from typing import Dict, Optional, Tuple, List

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

logger = logging.getLogger(__name__)


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
        effective_size_threshold : float, optional
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
        effective_size_threshold = kwargs.get('effective_size_threshold', None)
        
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
            effective_size_threshold=effective_size_threshold
        )
        
        # Save M matrices to each species
        print("\nSaving M matrices to species objects...")
        for species, M_matrix in results.items():
            self._species_data[species].M = M_matrix
            print(f"✓ Saved M matrix for {species}: {M_matrix.shape}")
        
        print("\nMulti-view ICA completed successfully!")
    
    def run_robust_multiview_ica(
        self, 
        a: Dict[str, int],
        c: int,
        num_runs: int = 100,
        mode: str = 'gpu',
        seed: int = 42,
        effective_size_threshold: Optional[float] = 5,
        effective_size_threshold_core: Optional[float] = None,
        effective_size_threshold_unique: Optional[float] = None,
        num_top_gene: int = 20,
        save_plots: Optional[str] = None
    ) -> Dict[str, pd.DataFrame]:
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
        effective_size_threshold : float, optional, default=5
            Cohen's d threshold for both core and unique components.
            Only used if specific thresholds are not provided.
        effective_size_threshold_core : float, optional
            Cohen's d threshold specifically for core components.
            If provided, overrides effective_size_threshold for core components.
        effective_size_threshold_unique : float, optional
            Cohen's d threshold specifically for unique components.
            If provided, overrides effective_size_threshold for unique components.
        num_top_gene : int, default=20
            Number of top genes to use when calculating Cohen's d effect size
        save_plots : str, optional
            Directory to save clustering plots. If None, no plots are saved.
            
        Returns
        -------
        dict
            Dictionary mapping species names to M matrices containing robust components.
            Each M matrix has columns for core components followed by unique components.
            
        Examples
        --------
        >>> # Run with optimized parameters
        >>> a_values = {'strain1': 50, 'strain2': 55, 'strain3': 60}
        >>> robust_results = multiModulon.run_robust_multiview_ica(
        ...     a=a_values,
        ...     c=30,
        ...     num_runs=100,
        ...     effective_size_threshold=5
        ... )
        
        >>> # Use different thresholds for core and unique components
        >>> robust_results = multiModulon.run_robust_multiview_ica(
        ...     a=a_values,
        ...     c=30,
        ...     effective_size_threshold_core=5,
        ...     effective_size_threshold_unique=3
        ... )
        """
        # Set effect size thresholds
        if effective_size_threshold_core is None:
            effective_size_threshold_core = effective_size_threshold if effective_size_threshold is not None else 5
        if effective_size_threshold_unique is None:
            effective_size_threshold_unique = effective_size_threshold if effective_size_threshold is not None else 5
            
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
        print(f"Effect size thresholds - Core: {effective_size_threshold_core}, Unique: {effective_size_threshold_unique}")
        
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
                    if effect_size >= effective_size_threshold_core:
                        # Apply sign convention
                        component = self._enforce_sign_convention(weight_vector)
                        core_components[species].append(component)
                
                # Process unique components  
                for comp_idx in range(c, a[species]):
                    weight_vector = M.iloc[:, comp_idx].values
                    
                    # Calculate effect size
                    effect_size = calculate_cohens_d_effect_size(weight_vector, seed, num_top_gene)
                    
                    # Only keep components above threshold
                    if effect_size >= effective_size_threshold_unique:
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
        
        # Print summary
        print(f"\n{'='*60}")
        print("Robust multi-view ICA completed!")
        print(f"{'='*60}")
        print(f"Total core clusters identified: {len(core_results['clusters'])}")
        for species in species_list:
            n_unique = len(unique_results[species])
            print(f"{species}: {n_unique} unique components")
        print(f"{'='*60}")
        
        return final_M_matrices
    
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