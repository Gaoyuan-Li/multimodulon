"""Core MultiModulon class for multi-species expression analysis."""

from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import json
import pickle
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, Optional, Tuple

from .utils.gff_parser import parse_gff
from .utils.bbh import BBHAnalyzer
from .utils.fasta_utils import extract_protein_sequences
from .multiview_ica import run_multiview_ica
from .multiview_ica_optimization import run_nre_optimization, calculate_cohens_d_effect_size

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
        self._log_tpm_norm = None
        self._X = None
        self._M = None
        self._A = None
        self._sample_sheet = None
        self._gene_table = None
    
    @property
    def log_tpm(self) -> pd.DataFrame:
        """Get log TPM expression matrix."""
        if self._log_tpm is None:
            self._load_log_tpm()
        return self._log_tpm
    
    @property
    def log_tpm_norm(self) -> pd.DataFrame:
        """Get normalized log TPM expression matrix."""
        if self._log_tpm_norm is None:
            self._load_log_tpm_norm()
        return self._log_tpm_norm
    
    @property
    def X(self) -> pd.DataFrame:
        """Get aligned expression matrix. Defaults to log_tpm_norm if not set."""
        if self._X is None:
            return self.log_tpm_norm
        return self._X
    
    @property
    def M(self) -> pd.DataFrame:
        """Get ICA mixing matrix (M matrix)."""
        if self._M is None:
            raise AttributeError(f"M matrix not yet computed for {self.species_name}. Run multiview ICA first.")
        return self._M
    
    @M.setter
    def M(self, value: pd.DataFrame):
        """Set ICA mixing matrix."""
        self._M = value
    
    @property
    def A(self) -> pd.DataFrame:
        """Get ICA activity matrix (A matrix)."""
        if self._A is None:
            raise AttributeError(f"A matrix not yet computed for {self.species_name}. Run generate_A() first.")
        return self._A
    
    @A.setter
    def A(self, value: pd.DataFrame):
        """Set ICA activity matrix."""
        self._A = value
    
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
        file_path = self.data_path / "expression_matrices" / "log_tpm.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"log_tpm.csv not found for {self.species_name}")
        
        # Read the file with first column as index
        self._log_tpm = pd.read_csv(file_path, index_col=0)
        logger.info(f"Loaded log_tpm for {self.species_name}: {self._log_tpm.shape}")
    
    def _load_log_tpm_norm(self):
        """Load normalized log TPM expression matrix."""
        file_path = self.data_path / "expression_matrices" / "log_tpm_norm.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"log_tpm_norm.csv not found for {self.species_name}")
        
        # Read the file with first column as index
        self._log_tpm_norm = pd.read_csv(file_path, index_col=0)
        logger.info(f"Loaded log_tpm_norm for {self.species_name}: {self._log_tpm_norm.shape}")
    
    def _load_sample_sheet(self):
        """Load sample metadata."""
        file_path = self.data_path / "samplesheet" / "samplesheet.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"samplesheet.csv not found for {self.species_name}")
        
        # Read the file with first column as index
        self._sample_sheet = pd.read_csv(file_path, index_col=0)
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
        # Check sample consistency between expression matrices and sample sheet
        log_tpm_samples = self.log_tpm.shape[1]
        log_tpm_norm_samples = self.log_tpm_norm.shape[1]
        sample_sheet_samples = self.sample_sheet.shape[0]
        
        # Check if number of samples match
        if log_tpm_samples != sample_sheet_samples:
            logger.warning(f"Sample count mismatch for {self.species_name}: "
                         f"log_tpm has {log_tpm_samples} samples, "
                         f"sample_sheet has {sample_sheet_samples} samples")
            return False
        
        if log_tpm_norm_samples != sample_sheet_samples:
            logger.warning(f"Sample count mismatch for {self.species_name}: "
                         f"log_tpm_norm has {log_tpm_norm_samples} samples, "
                         f"sample_sheet has {sample_sheet_samples} samples")
            return False
        
        if log_tpm_samples != log_tpm_norm_samples:
            logger.warning(f"Expression matrix mismatch for {self.species_name}: "
                         f"log_tpm has {log_tpm_samples} samples, "
                         f"log_tpm_norm has {log_tpm_norm_samples} samples")
            return False
        
        logger.info(f"Data validation passed for {self.species_name}")
        return True


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
    
    def generate_BBH(self, output_path: str = "Output_BBH", threads: int = 1):
        """
        Generate BBH files using existing protein.faa files from each strain.
        
        Parameters
        ----------
        output_path : str
            Path to save BBH results (default: "Output_BBH")
        threads : int
            Number of threads to use for BLAST (default: 1)
        """
        
        output_dir = Path(output_path)
        output_dir.mkdir(exist_ok=True)
        
        species_list = list(self._species_data.keys())
        logger.info(f"Generating BBH for {len(species_list)} species: {species_list}")
        
        # Use existing protein.faa files and create mapping from protein_id to locus_tag
        species_fasta_paths = {}
        protein_to_locus_maps = {}
        
        for species_name in species_list:
            species_data = self._species_data[species_name]
            ref_genome_dir = species_data.data_path / "ref_genome"
            
            # Check for existing protein.faa file
            protein_fasta = ref_genome_dir / "protein.faa"
            if not protein_fasta.exists():
                logger.error(f"Missing protein.faa file for {species_name} at {protein_fasta}")
                continue
            
            # Check for GFF file to create protein_id to locus_tag mapping
            gff_file = ref_genome_dir / "genomic.gff"
            if not gff_file.exists():
                logger.error(f"Missing genomic.gff file for {species_name}")
                continue
            
            # Create mapping from NCBI_GP (protein_id) to locus_tag
            protein_to_locus = self._create_protein_to_locus_mapping(gff_file, species_name)
            protein_to_locus_maps[species_name] = protein_to_locus
            species_fasta_paths[species_name] = protein_fasta
            
            logger.info(f"Found {len(protein_to_locus)} protein-to-locus mappings for {species_name}")
        
        # Check if BLAST tools are available
        try:
            subprocess.run(["makeblastdb", "-version"], capture_output=True, check=True)
            subprocess.run(["blastp", "-version"], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error("BLAST tools (makeblastdb, blastp) not found. Please install NCBI BLAST+:")
            logger.error("  Ubuntu/Debian: sudo apt-get install ncbi-blast+")
            logger.error("  macOS: brew install blast")
            logger.error("  Conda: conda install -c bioconda blast")
            raise RuntimeError("BLAST tools not installed")
        
        # First, run all BLAST comparisons to collect results
        blast_results = {}
        for species1 in species_list:
            for species2 in species_list:
                # Skip if we don't have protein files for these species
                if species1 not in species_fasta_paths or species2 not in species_fasta_paths:
                    logger.warning(f"Skipping {species1} vs {species2}: missing protein files")
                    continue
                
                key = (species1, species2)
                
                # Check if output file already exists
                output_file = output_dir / f"{species1}_vs_{species2}.csv"
                if output_file.exists():
                    logger.info(f"BBH file already exists: {output_file}")
                    # Load existing results
                    try:
                        blast_results[key] = pd.read_csv(output_file)
                    except Exception as e:
                        logger.warning(f"Failed to load existing file {output_file}: {e}")
                    continue
                
                logger.info(f"Running BLAST: {species1} vs {species2}")
                
                try:
                    # Run BLAST comparison
                    blast_df = self._run_blast_comparison_with_locus(
                        species_fasta_paths[species1],
                        species_fasta_paths[species2],
                        species1, species2,
                        protein_to_locus_maps[species1],
                        protein_to_locus_maps[species2],
                        threads
                    )
                    blast_results[key] = blast_df
                    
                except Exception as e:
                    logger.error(f"Failed to run BLAST for {species1} vs {species2}: {str(e)}")
                    blast_results[key] = pd.DataFrame(columns=['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 
                                                              'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 
                                                              'subjectEnd', 'eVal', 'bitScore', 'gene_length', 'COV', 'BBH'])
        
        # Now determine BBH for all pairs
        for species1 in species_list:
            for species2 in species_list:
                key = (species1, species2)
                if key not in blast_results:
                    continue
                
                output_file = output_dir / f"{species1}_vs_{species2}.csv"
                if output_file.exists():
                    continue
                
                df = blast_results[key].copy()
                
                if species1 == species2:
                    # Self-comparison, all hits are BBH
                    df['BBH'] = '<=>'
                else:
                    # Cross-species comparison, determine BBH
                    reverse_key = (species2, species1)
                    if reverse_key in blast_results:
                        reverse_df = blast_results[reverse_key]
                        # Create dictionaries for fast lookup
                        forward_hits = dict(zip(df['gene'], df['subject']))
                        reverse_hits = dict(zip(reverse_df['gene'], reverse_df['subject']))
                        
                        # Mark BBH
                        df['BBH'] = df.apply(
                            lambda row: '<=>' if (row['subject'] in reverse_hits and 
                                                reverse_hits[row['subject']] == row['gene']) else '',
                            axis=1
                        )
                    else:
                        logger.warning(f"No reverse comparison found for {species2} vs {species1}")
                        df['BBH'] = ''
                
                # Save results with index column
                df.reset_index(drop=True, inplace=True)
                df.to_csv(output_file, index=True, index_label='')
                logger.info(f"Saved BBH results to {output_file}")
        
        logger.info(f"BBH generation completed. Results saved to {output_dir}")
    
    def _run_blast_comparison(self, fasta1: Path, fasta2: Path, sp1: str, sp2: str) -> pd.DataFrame:
        """
        Run BLAST comparison between two species and return parsed results.
        
        Returns DataFrame with columns: gene, subject
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Copy FASTA files to temp directory
            tmp_fasta1 = tmpdir / f"{sp1}.faa"
            tmp_fasta2 = tmpdir / f"{sp2}.faa"
            
            subprocess.run(f"cp {fasta1} {tmp_fasta1}", shell=True, check=True)
            subprocess.run(f"cp {fasta2} {tmp_fasta2}", shell=True, check=True)
            
            # Create BLAST database
            subprocess.run([
                "makeblastdb", "-in", str(tmp_fasta2), "-dbtype", "prot", "-out", str(tmpdir / f"{sp2}_db")
            ], check=True, capture_output=True)
            
            # Run BLAST
            blast_output = tmpdir / "blast_results.txt"
            subprocess.run([
                "blastp", "-query", str(tmp_fasta1), "-db", str(tmpdir / f"{sp2}_db"),
                "-out", str(blast_output), "-outfmt", "6 qseqid sseqid evalue bitscore",
                "-max_target_seqs", "1", "-evalue", "1e-5"
            ], check=True, capture_output=True)
            
            # Parse BLAST results
            if blast_output.exists() and blast_output.stat().st_size > 0:
                df = pd.read_csv(
                    blast_output, sep='\\t',
                    names=['gene', 'subject', 'evalue', 'bitscore']
                )
                # Keep only gene and subject columns for BBH format
                return df[['gene', 'subject']]
            else:
                return pd.DataFrame(columns=['gene', 'subject'])
    
    def _create_protein_to_locus_mapping(self, gff_file: Path, species_name: str) -> dict:
        """
        Create mapping from protein_id (NCBI_GP) to locus_tag from GFF file.
        Handles different GFF formats including W3110 which uses Note field instead of locus_tag.
        
        For W3110 strain, extracts JW-numbers (e.g., JW4367) from Note field format "ECK0001:JW4367:b0001".
        Issues warnings when fallback extraction methods are used.
        
        Parameters
        ----------
        gff_file : Path
            Path to the GFF file
        species_name : str
            Name of the species/strain (folder name from Input_Data)
            
        Returns
        -------
        dict
            Mapping from protein_id to locus_tag
        """
        protein_to_locus = {}
        
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature_type = parts[2]
                if feature_type == 'CDS':
                    attributes = parts[8]
                    
                    # Parse attributes
                    attr_dict = {}
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attr_dict[key] = value
                    
                    # Get protein_id from NCBI_GP in Dbxref
                    protein_id = None
                    if 'Dbxref' in attr_dict:
                        dbxrefs = attr_dict['Dbxref'].split(',')
                        for dbxref in dbxrefs:
                            if dbxref.startswith('NCBI_GP:'):
                                protein_id = dbxref.replace('NCBI_GP:', '')
                                break
                    
                    # Get locus_tag - handle different formats
                    locus_tag = attr_dict.get('locus_tag', None)
                    
                    # If no locus_tag found, try to extract from Note field (W3110 format)
                    if not locus_tag and 'Note' in attr_dict:
                        # Issue warning for locus_tag extraction
                        warning_key = f"{species_name}_note_extraction"
                        if warning_key not in getattr(self, '_warned_extractions', set()):
                            print(f"WARNING: '{species_name}' - No locus_tag found in GFF file, extracting from Note field")
                            if not hasattr(self, '_warned_extractions'):
                                self._warned_extractions = set()
                            self._warned_extractions.add(warning_key)
                        
                        note_value = attr_dict['Note']
                        # Parse "ECK0001:JW4367:b0001" format - use the JW-number for W3110
                        if ':' in note_value:
                            parts = note_value.split(':')
                            # Look for JW-number first (preferred for W3110)
                            for part in parts:
                                if part.startswith('JW') and len(part) > 2:
                                    locus_tag = part
                                    break
                            # If no JW-number found, fall back to b-number
                            if not locus_tag:
                                for part in reversed(parts):
                                    if part.startswith('b') and len(part) > 1:
                                        locus_tag = part
                                        break
                    
                    # If still no locus_tag, try using gene name as fallback
                    if not locus_tag and 'gene' in attr_dict:
                        warning_key = f"{species_name}_gene_fallback"
                        if warning_key not in getattr(self, '_warned_extractions', set()):
                            print(f"WARNING: Species '{species_name}' - No locus_tag or Note field found in GFF file, using gene name as fallback")
                            if not hasattr(self, '_warned_extractions'):
                                self._warned_extractions = set()
                            self._warned_extractions.add(warning_key)
                        locus_tag = attr_dict['gene']
                    
                    # If still no locus_tag, use protein_id itself as last resort
                    if not locus_tag and protein_id:
                        warning_key = f"{species_name}_protein_id_fallback"
                        if warning_key not in getattr(self, '_warned_extractions', set()):
                            print(f"WARNING: Species '{species_name}' - No identifiable locus_tag in GFF file, using protein_id as last resort")
                            if not hasattr(self, '_warned_extractions'):
                                self._warned_extractions = set()
                            self._warned_extractions.add(warning_key)
                        locus_tag = protein_id
                    
                    if protein_id and locus_tag:
                        protein_to_locus[protein_id] = locus_tag
        
        return protein_to_locus
    
    def _run_blast_comparison_with_locus(self, fasta1: Path, fasta2: Path, sp1: str, sp2: str,
                                        protein_to_locus1: dict, protein_to_locus2: dict, threads: int = 1) -> pd.DataFrame:
        """
        Run BLAST comparison between two species and return results with locus_tags.
        
        Returns DataFrame with columns matching the expected BBH output format.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Copy FASTA files to temp directory
            tmp_fasta1 = tmpdir / f"{sp1}.faa"
            tmp_fasta2 = tmpdir / f"{sp2}.faa"
            
            subprocess.run(f"cp {fasta1} {tmp_fasta1}", shell=True, check=True)
            subprocess.run(f"cp {fasta2} {tmp_fasta2}", shell=True, check=True)
            
            # Create BLAST database
            subprocess.run([
                "makeblastdb", "-in", str(tmp_fasta2), "-dbtype", "prot", "-out", str(tmpdir / f"{sp2}_db")
            ], check=True, capture_output=True)
            
            # Run BLAST with extended output format
            blast_output = tmpdir / "blast_results.txt"
            subprocess.run([
                "blastp", "-query", str(tmp_fasta1), "-db", str(tmpdir / f"{sp2}_db"),
                "-out", str(blast_output), 
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen",
                "-max_target_seqs", "1", "-evalue", "1e-5",
                "-num_threads", str(threads)
            ], check=True, capture_output=True)
            
            # Parse BLAST results
            if blast_output.exists() and blast_output.stat().st_size > 0:
                df = pd.read_csv(
                    blast_output, sep='\t',
                    names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                           'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen']
                )
                
                # Convert protein IDs to locus_tags
                df['gene'] = df['qseqid'].map(protein_to_locus1)
                df['subject'] = df['sseqid'].map(protein_to_locus2)
                
                # Filter out rows where mapping failed
                df = df.dropna(subset=['gene', 'subject'])
                
                # Rename columns to match expected format
                df = df.rename(columns={
                    'pident': 'PID',
                    'length': 'alnLength',
                    'mismatch': 'mismatchCount',
                    'gapopen': 'gapOpenCount',
                    'qstart': 'queryStart',
                    'qend': 'queryEnd',
                    'sstart': 'subjectStart',
                    'send': 'subjectEnd',
                    'evalue': 'eVal',
                    'bitscore': 'bitScore',
                    'qlen': 'gene_length'
                })
                
                # Calculate coverage
                df['COV'] = df['alnLength'] / df['gene_length']
                
                # For self-comparisons, mark as bidirectional best hit
                if sp1 == sp2:
                    df['BBH'] = '<=>'  
                else:
                    # For cross-species comparisons, we need to run reciprocal BLAST to determine BBH
                    # For now, leave empty as the full BBH determination requires both directions
                    df['BBH'] = ''
                
                # Select and order columns to match expected format
                columns = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount',
                          'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore',
                          'gene_length', 'COV', 'BBH']
                
                return df[columns]
            else:
                return pd.DataFrame(columns=['gene', 'subject', 'PID', 'alnLength', 'mismatchCount',
                                           'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart',
                                           'subjectEnd', 'eVal', 'bitScore', 'gene_length', 'COV', 'BBH'])
    
    def align_genes(self, input_bbh_dir: str = "Output_BBH", output_dir: str = "Output_Gene_Info", 
                    reference_order: list[str] = None, bbh_threshold: float = None):
        """
        Align genes across all species using Union-Find algorithm to create combined gene database.
        
        Parameters
        ----------
        input_bbh_dir : str
            Path to the directory containing BBH CSV files (default: "Output_BBH")
        output_dir : str
            Path to save the combined gene database (default: "Output_Gene_Info")
        reference_order : list[str], optional
            List of strain names in the desired column order for the combined gene database.
            If not provided, strains will be ordered as they appear in the species data.
        bbh_threshold : float, optional
            Minimum percent identity (PID) threshold for BBH relationships.
            Gene pairs with PID below this threshold will be treated as different genes.
            If not provided, all BBH relationships are used regardless of PID.
        """
        from collections import defaultdict
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        bbh_dir = Path(input_bbh_dir)
        if not bbh_dir.exists():
            raise FileNotFoundError(f"BBH directory not found: {input_bbh_dir}")
        
        # Get list of strains
        if reference_order:
            # Validate that all strains in reference_order exist in species_data
            available_strains = set(self._species_data.keys())
            for strain in reference_order:
                if strain not in available_strains:
                    raise ValueError(f"Strain '{strain}' in reference_order not found in loaded species data. Available strains: {list(available_strains)}")
            strains = reference_order
        else:
            strains = list(self._species_data.keys())
        logger.info(f"Aligning genes for strains: {strains}")
        
        # Step 1: Collect all genes from each strain's self-comparison file
        strain_genes = {strain: set() for strain in strains}
        
        for strain in strains:
            self_file = bbh_dir / f'{strain}_vs_{strain}.csv'
            if self_file.exists():
                df = pd.read_csv(self_file)
                if 'gene' in df.columns:
                    genes = df['gene'].unique()
                    strain_genes[strain].update(genes)
                    logger.info(f"Found {len(genes)} genes for {strain}")
            else:
                # If self-comparison file doesn't exist, get genes from expression data
                species_data = self._species_data[strain]
                genes = species_data.log_tpm.index.tolist()
                strain_genes[strain].update(genes)
                logger.info(f"Using {len(genes)} genes from expression data for {strain}")
        
        # Step 2: Initialize Union-Find structure
        class UnionFind:
            def __init__(self):
                self.parent = {}
            
            def find(self, node):
                if node not in self.parent:
                    self.parent[node] = node
                while self.parent[node] != node:
                    self.parent[node] = self.parent[self.parent[node]]  # Path compression
                    node = self.parent[node]
                return node
            
            def union(self, node1, node2):
                root1 = self.find(node1)
                root2 = self.find(node2)
                if root1 != root2:
                    self.parent[root2] = root1
        
        uf = UnionFind()
        
        # Step 3: Process all pairwise BBH files to build equivalence relationships
        if bbh_dir.exists():
            for filename in bbh_dir.iterdir():
                if not filename.name.endswith('.csv'):
                    continue
                # Parse strain names from filename
                parts = filename.stem.split('_vs_')
                if len(parts) != 2:
                    continue
                x, y = parts[0], parts[1]
                
                if x not in strains or y not in strains:
                    continue
                
                try:
                    df = pd.read_csv(filename)
                    if 'gene' in df.columns and 'subject' in df.columns and 'BBH' in df.columns:
                        # Only process rows with BBH marker '<=>'
                        bbh_rows = df[df['BBH'] == '<=>']
                        logger.info(f"Processing {len(bbh_rows)} BBH relationships from {filename.name}")
                        
                        for _, row in bbh_rows.iterrows():
                            gene_x = (x, row['gene'])
                            gene_y = (y, row['subject'])
                            # Ensure genes exist in their respective strain's self-reported list
                            if row['gene'] in strain_genes[x] and row['subject'] in strain_genes[y]:
                                # Check PID threshold if provided
                                if bbh_threshold is not None and 'PID' in df.columns:
                                    if row['PID'] >= bbh_threshold:
                                        uf.union(gene_x, gene_y)
                                    # If PID is below threshold, genes are treated as different (no union)
                                else:
                                    # No threshold specified or PID column missing, use all relationships
                                    uf.union(gene_x, gene_y)
                    else:
                        logger.warning(f"Required columns missing in {filename}: {df.columns.tolist()}")
                except Exception as e:
                    logger.warning(f"Failed to process {filename}: {str(e)}")
        
        # Step 4: Group all genes into their connected components
        components = defaultdict(lambda: {strain: [] for strain in strains})
        for strain in strains:
            for gene in strain_genes[strain]:
                node = (strain, gene)
                root = uf.find(node)
                components[root][strain].append(gene)
        
        # Step 5: Create the final DataFrame
        rows = []
        for comp in components.values():
            row = {}
            for strain in strains:
                genes = comp[strain]
                if len(genes) == 1:
                    row[strain] = genes[0]
                elif len(genes) > 1:
                    # Handle unexpected duplicates: take the first occurrence
                    row[strain] = genes[0]
                else:
                    row[strain] = None
            rows.append(row)
        
        # Create DataFrame and sort by leftmost non-None value
        df = pd.DataFrame(rows)
        
        # Add sorting logic: prioritize rows with non-None values in leftmost columns
        def get_sort_key(row):
            # Create a tuple for sorting: (column_index_of_first_non_none, value_at_that_position)
            for i, strain in enumerate(strains):
                val = row[strain]
                if pd.notna(val) and val != "None" and val is not None:
                    return (i, val)
            return (len(strains), "")  # Rows with all None values go to the end
        
        df['_sort_key'] = df.apply(get_sort_key, axis=1)
        
        # Sort by column index first, then by value within each column group
        df = df.sort_values('_sort_key').drop('_sort_key', axis=1).reset_index(drop=True)
        
        output_file = output_path / "combined_gene_db.csv"
        df.to_csv(output_file, index=False)
        
        # Store the combined gene database
        self._combined_gene_db = df
        
        logger.info(f"Gene alignment completed. Combined gene database saved to {output_file}")
        logger.info(f"Total gene groups: {len(df)}")
        
        # Create aligned expression matrices
        self._create_aligned_expression_matrices(df, strains)
        
        return df
    
    def _create_aligned_expression_matrices(self, combined_gene_db: pd.DataFrame, strains: list[str]):
        """
        Create aligned expression matrices (X) for all strains with same row indexes.
        """
        logger.info("Creating aligned expression matrices...")
        
        # Define a helper function to find the first valid gene ID for row labeling
        def find_first_valid_gene(row):
            for strain in strains:
                val = row[strain]
                # Treat "None" strings, actual NaN, or None as invalid
                if pd.notna(val) and val != "None" and val is not None:
                    return val
            return None
        
        # Apply this helper to each row of combined_gene_db
        combined_gene_db["row_label"] = combined_gene_db.apply(find_first_valid_gene, axis=1)
        
        # Remove rows with no valid genes
        valid_rows = combined_gene_db.dropna(subset=['row_label'])
        
        # The combined_gene_db is already sorted by leftmost non-None value, so use it as is
        df_sorted = valid_rows.copy()
        
        # Set the gene IDs as the new index
        df_sorted = df_sorted.set_index("row_label")
        
        # Get original expression data for each strain
        expression_dfs = {}
        for strain in strains:
            species_data = self._species_data[strain]
            expression_dfs[strain] = species_data.log_tpm_norm  # Use the normalized expression data
        
        # Build aligned expression matrices
        for strain in strains:
            # Prepare an empty DataFrame with same columns as original but reindexed
            if strain in expression_dfs and not expression_dfs[strain].empty:
                cols = expression_dfs[strain].columns
                new_df = pd.DataFrame(index=df_sorted.index, columns=cols, dtype=float)
                
                # For each row in df_sorted, find the gene ID that corresponds to this strain
                for gene_label in new_df.index:
                    strain_gene_id = df_sorted.loc[gene_label, strain]
                    if pd.notna(strain_gene_id) and strain_gene_id != "None" and strain_gene_id is not None:
                        # If the gene is in the strain's expression file, copy it. If not, fill with 0.
                        if strain_gene_id in expression_dfs[strain].index:
                            new_df.loc[gene_label] = expression_dfs[strain].loc[strain_gene_id]
                        else:
                            new_df.loc[gene_label] = 0
                    else:
                        new_df.loc[gene_label] = 0
                
                # Update only the X matrix (preserve original log_tpm_norm)
                self._species_data[strain]._X = new_df
                logger.info(f"Created aligned expression matrix for {strain}: {new_df.shape}")
            else:
                logger.warning(f"No expression data found for {strain}")
    
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
        **kwargs : dict
            Arguments for multi-view ICA:
            - a1, a2, ..., a6 : int
                Number of independent components for each species/strain
            - c : int
                Number of core components across all species
            - mode : str, optional
                'gpu' or 'cpu' mode (default: 'gpu')
            - effective_size_threshold : float, optional
                Cohen's d effect size threshold for filtering components (comparing top 20 genes vs rest).
                If provided, only keeps components above this threshold.
        
        Returns
        -------
        None
            Updates each species' M matrix in place
        
        Examples
        --------
        # For 6 species (traditional way)
        >>> multiModulon.run_multiview_ica(a1=50, a2=50, a3=50, a4=50, a5=50, a6=50, c=30)
        
        # For any number of species (using list)
        >>> multiModulon.run_multiview_ica(a=[50, 50, 50], c=30)  # for 3 species
        
        # For any number of species (same a for all)
        >>> multiModulon.run_multiview_ica(a=50, c=30)
        
        # With Cohen's d threshold filtering (top 20 genes vs rest)
        >>> multiModulon.run_multiview_ica(a=50, c=30, effective_size_threshold=0.2)
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
    
    def optimize_number_of_core_components(
        self,
        max_k: Optional[int] = None,
        step: int = 5,
        max_a_per_view: Optional[int] = None,
        train_frac: float = 0.75,
        num_runs: int = 1,
        mode: str = 'gpu',
        seed: int = 42,
        save_plot: Optional[str] = None,
        metric: str = 'nre',
        threshold: Optional[float] = None,
        effective_size_threshold: float = 5,
        num_top_gene: int = 20
    ) -> Tuple[int, Dict[int, float]]:
        """
        Optimize the number of core components using specified metric.
        
        This method uses cross-validation to find the optimal number of core
        components (k) by optimizing the selected metric.
        
        **Available Metrics:**
        
        1. **NRE (Normalized Reconstruction Error):**
           - Measures how much the "core" components differ across views
           - When k < k_true: Components are truly shared → low residuals → low NRE
           - When k = k_true: Optimal sharing → minimum NRE
           - When k > k_true: Including view-specific components → higher NRE
        
        2. **effect_size (Cohen's d Effect Size):**
           - Calculates Cohen's d between top 20 genes and remaining genes
           - Effect Size = |mean_top20 - mean_rest| / pooled_std
           - Measures separation between highly weighted genes and the rest
           - Higher effect sizes indicate clearer gene membership
           - Runs square ICA (a=c=k) to get mixing matrices M
        
        Parameters
        ----------
        max_k : int, optional
            Maximum k to test. If None, uses the maximum dimension from generate_X
        step : int, default=5
            Step size for k candidates (tests k = step, 2*step, 3*step, ...)
        max_a_per_view : int, optional
            Maximum components per view. If None, uses max_k
        train_frac : float, default=0.75
            Fraction of data to use for training
        num_runs : int, default=1
            Number of cross-validation runs
        mode : str, default='gpu'
            'gpu' or 'cpu' mode
        seed : int, default=42
            Random seed for reproducibility
        save_plot : str, optional
            Path to save the metric vs k plot. If None, displays the plot
        metric : str, default='nre'
            Optimization metric: 'nre' or 'effect_size' (Cohen's d effect size)
        threshold : float, optional
            For effect_size metric: Cohen's d threshold for effect size analysis and visualization
        effective_size_threshold : float, default=5
            Minimum Cohen's d effect size threshold for counting components.
            Only components with effect size above this threshold are counted as "effective".
        num_top_gene : int, default=20
            Number of top genes to use when calculating Cohen's d effect size
        
        Returns
        -------
        best_k : int
            Optimal number of core components
        metric_scores : dict
            Dictionary mapping k values to mean metric scores (or number of above-threshold components for effect_size metric)
        
        Examples
        --------
        >>> # Using NRE metric (default)
        >>> best_k, nre_scores = multiModulon.optimize_number_of_core_components()
        >>> print(f"Optimal number of core components: {best_k}")
        
        >>> # Using Cohen's d effect size metric for regulatory structure analysis
        >>> best_k, effect_size_scores = multiModulon.optimize_number_of_core_components(metric='effect_size')
        >>> print(f"Optimal number of core components: {best_k}")
        
        >>> # Using Cohen's d with threshold analysis
        >>> best_k, effect_size_scores = multiModulon.optimize_number_of_core_components(metric='effect_size', threshold=0.1)
        >>> print(f"Optimal number of core components: {best_k}")
        """
        # Check prerequisites
        species_list = list(self._species_data.keys())
        n_species = len(species_list)
        if n_species < 2:
            raise ValueError(
                f"NRE optimization requires at least 2 species/strains, found {n_species}"
            )
        
        print(f"Optimizing core components for {n_species} species/strains: {species_list}")
        
        # Check if all species have X matrices
        for species in species_list:
            if self._species_data[species]._X is None:
                raise ValueError(
                    f"X matrix not found for {species}. "
                    "Please run generate_X() first to create aligned expression matrices."
                )
        
        # Determine max_k if not provided
        if max_k is None:
            # Get number of samples for each species (should be consistent after alignment)
            n_samples = []
            for species in species_list:
                X = self._species_data[species].X
                n_samples.append(X.shape[1])  # number of samples (columns)
            min_samples = min(n_samples)
            # Set max_k to the largest multiple of step that's less than min_samples
            # This ensures k candidates don't exceed data constraints
            max_k = ((min_samples - 1) // step) * step
            max_k = max(step, max_k)  # Ensure at least one k candidate
            print(f"Auto-determined max_k = {max_k} based on minimum samples ({min_samples})")
        
        # Generate k candidates
        k_candidates = list(range(step, min(max_k + 1, 100), step))
        if not k_candidates:
            raise ValueError(f"No valid k candidates with step={step} and max_k={max_k}")
        
        # Set max_a_per_view if not provided
        if max_a_per_view is None:
            max_a_per_view = max_k
        
        # Prepare X matrices
        species_X_matrices = {}
        for species in species_list:
            species_X_matrices[species] = self._species_data[species].X
        
        # Run optimization
        best_k, metric_scores, all_metric_per_k, fig = run_nre_optimization(
            species_X_matrices=species_X_matrices,
            k_candidates=k_candidates,
            max_a_per_view=max_a_per_view,
            train_frac=train_frac,
            num_runs=num_runs,
            mode=mode,
            seed=seed,
            metric=metric,
            threshold=threshold,
            effective_size_threshold=effective_size_threshold,
            num_top_gene=num_top_gene
        )
        
        # Save or display plot
        if save_plot:
            fig.savefig(save_plot, dpi=300, bbox_inches='tight')
            print(f"\nPlot saved to: {save_plot}")
        else:
            import matplotlib.pyplot as plt
            plt.show()
        
        # Store the optimal k for later use
        self._optimal_k = best_k
        
        return best_k, metric_scores
    
    def _check_component_consistency(self, A: pd.DataFrame, sample_sheet: pd.DataFrame, component_name: str) -> bool:
        """
        Check if a component is consistent across biological replicates.
        
        Returns True if component is consistent, False if it should be dropped.
        """
        # Build replicate groups
        if "biological_replicate" not in sample_sheet.columns:
            # If no biological_replicate column, assume all samples are independent
            return True
        
        st = sample_sheet.copy()
        grp_id = 0
        replicate_groups = []
        for v in st["biological_replicate"]:
            if v == 1 and replicate_groups:
                grp_id += 1
            replicate_groups.append(grp_id)
        st["replicate_group"] = replicate_groups
        sample_to_grp = st["replicate_group"].to_dict()
        
        # Get component values
        row = A.loc[component_name]
        
        # Check each replicate group
        for gid in st["replicate_group"].unique():
            cols = [c for c in A.columns if c in sample_to_grp and sample_to_grp[c] == gid]
            if not cols:
                continue
                
            vals = row[cols].astype(float).values
            if vals.size == 0:
                continue
            
            abs_vals = np.abs(vals)
            mixed_signs = (vals > 0).any() and (vals < 0).any()
            all_big = (abs_vals > 5).all()
            wide_span = abs_vals.max() - abs_vals.min() > 20
            
            if (mixed_signs and all_big) or wide_span:
                return False  # Inconsistent
        
        return True  # Consistent
    
    def optimize_number_of_unique_components(
        self,
        best_k: Optional[int] = None,
        step: int = 5,
        mode: str = 'gpu',
        seed: int = 42,
        save_plots: Optional[str] = None,
        effective_size_threshold: float = 5,
        num_top_gene: int = 20
    ) -> Dict[str, int]:
        """
        Optimize the number of unique components for each species.
        
        This method runs ICA with varying numbers of unique components for each species
        while keeping other species at c=k_best, and finds the optimal number that
        maximizes consistent unique components.
        
        Parameters
        ----------
        best_k : int, optional
            Number of core components. If None, uses the value from optimize_number_of_core_components
        step : int, default=5
            Step size for testing unique components (tests k_best+5, k_best+10, ...)
        mode : str, default='gpu'
            'gpu' or 'cpu' mode
        seed : int, default=42
            Random seed for reproducibility
        save_plots : str, optional
            Directory to save plots. If None, displays the plots
        effective_size_threshold : float, default=5
            Minimum Cohen's d effect size threshold for counting components
        num_top_gene : int, default=20
            Number of top genes to use when calculating Cohen's d effect size
            
        Returns
        -------
        dict
            Dictionary mapping species names to optimal number of components (a values)
        """
        from tqdm import tqdm
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
        
        # Use stored optimal k if not provided
        if best_k is None:
            if hasattr(self, '_optimal_k'):
                best_k = self._optimal_k
            else:
                raise ValueError("best_k not provided and no optimal k found. Run optimize_number_of_core_components first.")
        
        print(f"\nOptimizing unique components with core k = {best_k}")
        
        species_list = list(self._species_data.keys())
        n_species = len(species_list)
        
        # Prepare X matrices
        species_X_matrices = {}
        for species in species_list:
            species_X_matrices[species] = self._species_data[species].X
        
        # Results storage
        optimal_a_values = {}
        
        # Process each species
        for target_species in species_list:
            print(f"\n{'='*60}")
            print(f"Optimizing unique components for {target_species}")
            print(f"{'='*60}")
            
            # Get minimum dimension (number of samples)
            min_dim = species_X_matrices[target_species].shape[1]
            
            # Generate a candidates: from k_best to closest multiple of 5 below min_dim
            a_candidates = []
            a = best_k
            while a < min_dim - 5:
                a_candidates.append(a)
                a += step
            # Add the last valid a
            if a_candidates[-1] < min_dim - 5:
                a_candidates.append(((min_dim - 5) // step) * step)
            
            if not a_candidates:
                a_candidates = [best_k]
            
            print(f"Testing a values: {a_candidates}")
            
            # Store results
            consistent_counts = {}
            
            # Test each a value
            for a_test in tqdm(a_candidates, desc=f"Testing a values for {target_species}"):
                # Set up a_values: a_test for target species, best_k for others
                a_values = {}
                for species in species_list:
                    if species == target_species:
                        a_values[species] = a_test
                    else:
                        a_values[species] = best_k
                
                # Run ICA
                M_matrices = run_multiview_ica(
                    species_X_matrices,
                    a_values,
                    best_k,  # c = best_k
                    mode=mode,
                    seed=seed
                )
                
                # Generate A matrix for target species
                M = M_matrices[target_species]
                X = species_X_matrices[target_species]
                A = M.T @ X
                A.index = M.columns
                A.columns = X.columns
                
                # Get sample sheet
                sample_sheet = self._species_data[target_species].sample_sheet
                
                # Count consistent unique components
                consistent_unique = 0
                
                for comp in A.index:
                    if comp.startswith("Unique"):
                        # Get component index (Unique components start after core components)
                        comp_idx = int(comp.split('_')[1]) - 1 + best_k
                        
                        # Check Cohen's d threshold
                        weight_vector = M.iloc[:, comp_idx].values
                        effect_size = calculate_cohens_d_effect_size(weight_vector, seed, num_top_gene)
                        
                        if effect_size >= effective_size_threshold:
                            # Check consistency
                            if self._check_component_consistency(A, sample_sheet, comp):
                                consistent_unique += 1
                
                consistent_counts[a_test] = consistent_unique
            
            # Find optimal a using 90% growth method
            a_sorted = sorted(consistent_counts.keys())
            counts = [consistent_counts[a] for a in a_sorted]
            
            if len(a_sorted) < 2:
                optimal_a = a_sorted[0]
            else:
                max_count = max(counts)
                min_count = min(counts)
                count_range = max_count - min_count
                
                if count_range == 0:
                    optimal_a = a_sorted[0]
                else:
                    # Find where we reach 90% of the range
                    threshold_percentage = 0.9
                    target_count = min_count + threshold_percentage * count_range
                    
                    # Find the smallest a that achieves at least the target count
                    optimal_a = a_sorted[0]
                    for i, (a, count) in enumerate(zip(a_sorted, counts)):
                        if count >= target_count:
                            optimal_a = a
                            break
            
            optimal_a_values[target_species] = optimal_a
            
            # Create plot
            fig, ax = plt.subplots(figsize=(10, 6))
            
            ax.plot(a_sorted, counts, 'bo-', linewidth=2, markersize=8, 
                   label='Consistent unique components')
            
            # Highlight optimal a
            ax.axvline(x=optimal_a, color='red', linestyle='--', alpha=0.7, linewidth=2,
                      label=f'Optimal a = {optimal_a}')
            ax.scatter([optimal_a], [consistent_counts[optimal_a]], color='red', s=120, zorder=5,
                      edgecolors='darkred', linewidth=2)
            
            ax.set_xlabel('Number of Total Components (a)', fontsize=12)
            ax.set_ylabel('Number of Consistent Unique Components', fontsize=12)
            ax.set_title(f'Optimization of Unique Components for {target_species}', fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Force integer y-axis
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            
            # Add text box with result
            label_text = f'Optimal a = {optimal_a}\n{consistent_counts[optimal_a]} consistent unique components'
            ax.text(0.02, 0.98, label_text, 
                   transform=ax.transAxes, fontsize=11, fontweight='bold',
                   bbox=dict(boxstyle="round,pad=0.4", facecolor='lightblue', alpha=0.8),
                   verticalalignment='top')
            
            plt.tight_layout()
            
            # Save or display plot
            if save_plots:
                Path(save_plots).mkdir(exist_ok=True)
                plot_path = Path(save_plots) / f"unique_optimization_{target_species}.png"
                fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                print(f"Plot saved to: {plot_path}")
            else:
                plt.show()
            
            print(f"\nOptimal a for {target_species}: {optimal_a} ({consistent_counts[optimal_a]} consistent unique components)")
        
        print(f"\n{'='*60}")
        print("Optimization Summary")
        print(f"{'='*60}")
        print(f"Core components (c): {best_k}")
        for species, a in optimal_a_values.items():
            print(f"{species}: a = {a} (unique components: {a - best_k})")
        
        return optimal_a_values
            
        if self._bbh is not None:
            print(f"\nBBH analysis completed:")
            print(f"  - Total species pairs: {len(self._bbh) // 2}")
            
            # Print some BBH statistics
            for (sp1, sp2), bbh_df in list(self._bbh.items())[:3]:
                if not bbh_df.empty:
                    print(f"  - {sp1} vs {sp2}: {len(bbh_df)} orthologs")