"""Core MultiModulon class for multi-species expression analysis."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional, Any, List, Tuple
import logging
import json
import pickle
from collections import defaultdict

from .utils.gff_parser import parse_gff
from .utils.bbh import BBHAnalyzer
from .utils.fasta_utils import extract_protein_sequences

logger = logging.getLogger(__name__)


class SpeciesData:
    """Container for single species data."""
    
    def __init__(self, species_name: str, data_path: Path):
        """
        Initialize species data container.
        
        Parameters
        ----------
        species_name : str
            Name of the species
        data_path : Path
            Path to the species data directory
        """
        self.species_name = species_name
        self.data_path = data_path
        self._log_tpm = None
        self._X = None
        self._sample_sheet = None
        self._gene_table = None
    
    @property
    def log_tpm(self) -> pd.DataFrame:
        """Get log TPM expression matrix."""
        if self._log_tpm is None:
            self._load_log_tpm()
        return self._log_tpm
    
    @property
    def X(self) -> pd.DataFrame:
        """Get normalized log TPM expression matrix."""
        if self._X is None:
            self._load_X()
        return self._X
    
    @property
    def sample_sheet(self) -> pd.DataFrame:
        """Get sample metadata."""
        if self._sample_sheet is None:
            self._load_sample_sheet()
        return self._sample_sheet
    
    @property
    def gene_table(self) -> pd.DataFrame:
        """Get gene annotation table."""
        if self._gene_table is None:
            self._load_gene_table()
        return self._gene_table
    
    def _load_log_tpm(self):
        """Load log TPM expression matrix."""
        file_path = self.data_path / "expression_matrices" / "log_tpm.tsv"
        if not file_path.exists():
            raise FileNotFoundError(f"log_tpm.tsv not found for {self.species_name}")
        
        # Read the file with first column as index
        self._log_tpm = pd.read_csv(file_path, sep='\t', index_col=0)
        logger.info(f"Loaded log_tpm for {self.species_name}: {self._log_tpm.shape}")
    
    def _load_X(self):
        """Load normalized log TPM expression matrix."""
        file_path = self.data_path / "expression_matrices" / "log_tpm_norm.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"log_tpm_norm.csv not found for {self.species_name}")
        
        # Read the file with first column as index
        self._X = pd.read_csv(file_path, index_col=0)
        logger.info(f"Loaded X (normalized) for {self.species_name}: {self._X.shape}")
    
    def _load_sample_sheet(self):
        """Load sample metadata."""
        file_path = self.data_path / "samplesheet" / "samplesheet.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"samplesheet.csv not found for {self.species_name}")
        
        # Read the file and set appropriate index
        self._sample_sheet = pd.read_csv(file_path)
        # Set the sample column as index if it exists
        if 'sample' in self._sample_sheet.columns:
            self._sample_sheet = self._sample_sheet.set_index('sample')
        logger.info(f"Loaded sample_sheet for {self.species_name}: {self._sample_sheet.shape}")
    
    def _load_gene_table(self):
        """Load gene annotation table from GFF file."""
        gff_file = self.data_path / "ref_genome" / "genomic.gff"
        if not gff_file.exists():
            raise FileNotFoundError(f"genomic.gff not found for {self.species_name}")
        
        self._gene_table = parse_gff(gff_file)
        logger.info(f"Loaded gene_table for {self.species_name}: {self._gene_table.shape}")
    
    def validate_data(self):
        """Validate data consistency."""
        # Check if gene IDs in log_tpm match gene_table
        log_tpm_genes = set(self.log_tpm.index)
        gene_table_genes = set(self.gene_table.index)
        
        # Find genes in log_tpm
        log_tpm_gene_prefixes = {g.replace('gene-', '') for g in log_tpm_genes if g.startswith('gene-')}
        
        # Check overlap
        common_genes = log_tpm_gene_prefixes.intersection(gene_table_genes)
        
        if len(common_genes) == 0:
            logger.warning(f"No common genes found between log_tpm and gene_table for {self.species_name}")
        else:
            logger.info(f"Found {len(common_genes)} common genes for {self.species_name}")
        
        return len(common_genes) > 0


class MultiModulon:
    """
    Main class for multi-species expression analysis.
    
    This class provides a unified interface for accessing expression data,
    gene annotations, and ortholog relationships across multiple species.
    """
    
    def __init__(self, data_folder: str, run_bbh: bool = True):
        """
        Initialize MultiModulon object.
        
        Parameters
        ----------
        data_folder : str
            Path to the folder containing species data
        run_bbh : bool
            Whether to run BBH analysis during initialization
        """
        self.data_folder = Path(data_folder)
        if not self.data_folder.exists():
            raise ValueError(f"Data folder not found: {data_folder}")
        
        # Initialize storage
        self._species_data: Dict[str, SpeciesData] = {}
        self._bbh: Optional[Dict[Tuple[str, str], pd.DataFrame]] = None
        
        # Load species data
        self._load_species_data()
        
        # Run BBH analysis if requested
        if run_bbh:
            self._run_bbh_analysis()
    
    def _load_species_data(self):
        """Load data for all species in the data folder."""
        # Find all species directories
        species_dirs = [d for d in self.data_folder.iterdir() 
                       if d.is_dir() and not d.name.startswith('.')]
        
        if not species_dirs:
            raise ValueError(f"No species directories found in {self.data_folder}")
        
        logger.info(f"Found {len(species_dirs)} species directories")
        
        # Load each species
        for species_dir in species_dirs:
            species_name = species_dir.name
            logger.info(f"Loading data for {species_name}")
            
            try:
                species_data = SpeciesData(species_name, species_dir)
                # Validate data
                if species_data.validate_data():
                    self._species_data[species_name] = species_data
                else:
                    logger.warning(f"Skipping {species_name} due to validation failure")
            except Exception as e:
                logger.error(f"Failed to load {species_name}: {str(e)}")
    
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
    def bbh(self) -> Optional[Dict[Tuple[str, str], pd.DataFrame]]:
        """Get BBH results."""
        return self._bbh
    
    @property
    def species(self) -> List[str]:
        """Get list of loaded species."""
        return list(self._species_data.keys())
    
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
    
    def summary(self):
        """Print summary of loaded data."""
        print(f"MultiModulon Summary")
        print(f"====================")
        print(f"Data folder: {self.data_folder}")
        print(f"Number of species: {len(self._species_data)}")
        print(f"\nSpecies loaded:")
        
        for species_name, species_data in self._species_data.items():
            print(f"\n  {species_name}:")
            print(f"    - Genes: {len(species_data.gene_table)}")
            print(f"    - Samples: {species_data.log_tpm.shape[1]}")
            print(f"    - Expression matrix shape: {species_data.log_tpm.shape}")
            
        if self._bbh is not None:
            print(f"\nBBH analysis completed:")
            print(f"  - Total species pairs: {len(self._bbh) // 2}")
            
            # Print some BBH statistics
            for (sp1, sp2), bbh_df in list(self._bbh.items())[:3]:
                if not bbh_df.empty:
                    print(f"  - {sp1} vs {sp2}: {len(bbh_df)} orthologs")