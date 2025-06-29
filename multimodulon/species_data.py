"""Species data container for MultiModulon analysis."""

from __future__ import annotations

import pandas as pd
from pathlib import Path
import logging

from .utils.gff_parser import parse_gff

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