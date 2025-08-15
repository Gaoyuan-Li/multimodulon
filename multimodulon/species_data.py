"""Species data container for MultiModulon analysis."""

from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path
import logging

from .gff_utils import gff2pandas

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
        self._M_thresholds = None
        self._presence_matrix = None
    
    @property
    def log_tpm(self) -> pd.DataFrame:
        """Get log TPM expression matrix."""
        if self._log_tpm is None:
            self._load_log_tpm()
        return self._log_tpm
    
    @log_tpm.setter
    def log_tpm(self, value: pd.DataFrame):
        """Set log TPM expression matrix."""
        self._log_tpm = value
    
    @property
    def log_tpm_norm(self) -> pd.DataFrame:
        """Get normalized log TPM expression matrix."""
        if self._log_tpm_norm is None:
            self._load_log_tpm_norm()
        return self._log_tpm_norm
    
    @log_tpm_norm.setter
    def log_tpm_norm(self, value: pd.DataFrame):
        """Set normalized log TPM expression matrix."""
        self._log_tpm_norm = value
    
    @property
    def X(self) -> pd.DataFrame:
        """Get aligned expression matrix. Defaults to log_tpm_norm if not set."""
        if self._X is None:
            return self.log_tpm_norm
        return self._X
    
    @X.setter
    def X(self, value: pd.DataFrame):
        """Set aligned expression matrix."""
        self._X = value
    
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
    
    @sample_sheet.setter
    def sample_sheet(self, value: pd.DataFrame):
        """Set sample metadata."""
        self._sample_sheet = value
    
    @property
    def gene_table(self) -> pd.DataFrame:
        """Get gene annotation table."""
        if self._gene_table is None:
            self._load_gene_table()
        return self._gene_table
    
    @gene_table.setter
    def gene_table(self, value: pd.DataFrame):
        """Set gene annotation table."""
        self._gene_table = value
    
    @property
    def M_thresholds(self) -> pd.DataFrame:
        """Get M matrix thresholds."""
        if self._M_thresholds is None:
            raise AttributeError(f"M_thresholds not yet computed for {self.species_name}. Run optimize_M_thresholds() first.")
        return self._M_thresholds
    
    @M_thresholds.setter
    def M_thresholds(self, value: pd.DataFrame):
        """Set M matrix thresholds."""
        self._M_thresholds = value
    
    @property
    def presence_matrix(self) -> pd.DataFrame:
        """Get presence matrix (binarized M matrix based on thresholds)."""
        if self._presence_matrix is None:
            raise AttributeError(f"presence_matrix not yet computed for {self.species_name}. Run optimize_M_thresholds() first.")
        return self._presence_matrix
    
    @presence_matrix.setter
    def presence_matrix(self, value: pd.DataFrame):
        """Set presence matrix."""
        self._presence_matrix = value
    
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
        
        self._gene_table = gff2pandas(str(gff_file), index='locus_tag')
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
    
    def calculate_explained_variance(
        self, genes=None, samples=None, imodulons=None
    ):
        """
        Calculate the explained variance for each component sorted by importance.
        
        This method computes the explained variance based on the reconstruction
        of the expression matrix from the product of M and A matrices, following
        the ICA decomposition X = MA.
        
        Parameters
        ----------
        genes : list, str, or None, optional
            Specific genes to include in the calculation. If None, all genes are used.
            Can be a single gene (str) or list of genes.
        samples : list, str, or None, optional
            Specific samples to include in the calculation. If None, all samples are used.
            Can be a single sample (str) or list of samples.
        imodulons : list, str, int, or None, optional
            Specific iModulons to include in the calculation. If None, all iModulons are used.
            Can be a single iModulon (str/int) or list of iModulons.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with columns 'iModulon' and 'Explained Variance',
            sorted by explained variance in descending order.
        
        Raises
        ------
        AttributeError
            If M or A matrices are not available.
        
        Examples
        --------
        >>> # Get explained variance for all components
        >>> exp_var = multimodulon['MG1655'].calculate_explained_variance()
        >>> print(exp_var.head())
           iModulon  Explained Variance
        0        5            0.052134
        1       12            0.041235
        2        3            0.038567
        """
        # Get the log_tpm matrix (use X if available, otherwise log_tpm_norm)
        log_tpm = self.X if self._X is not None else self.log_tpm_norm
        
        # Check inputs
        if genes is None:
            genes = log_tpm.index
        elif isinstance(genes, str):
            genes = [genes]

        if samples is None:
            samples = log_tpm.columns
        elif isinstance(samples, str):
            samples = [samples]

        if imodulons is None:
            imodulons = self.M.columns
        elif isinstance(imodulons, str) or isinstance(imodulons, int):
            imodulons = [imodulons]

        centered = log_tpm
        
        # Account for normalization procedures before ICA (X=SA-x_mean)
        baseline = centered.subtract(centered.mean(axis=0), axis=1)
        baseline = baseline.loc[genes, samples]

        # Initialize variables
        base_err = np.linalg.norm(baseline) ** 2
        MA = np.zeros(baseline.shape)
        rec_var = [0]
        ma_arrs = {}
        ma_weights = {}
        explained_variance_dict = {}
        i = 0
        
        # Get individual modulon contributions
        for k in imodulons:
            ma_arr = np.dot(
                self.M.loc[genes, k].values.reshape(len(genes), 1),
                self.A.loc[k, samples].values.reshape(1, len(samples)),
            )
            ma_arrs[k] = ma_arr
            ma_weights[k] = np.sum(ma_arr**2)

        # Sum components in order of most important component first
        sorted_mods = sorted(ma_weights, key=ma_weights.get, reverse=True)
        
        # Compute reconstructed variance
        for k in sorted_mods:
            MA = MA + ma_arrs[k]
            sa_err = np.linalg.norm(MA - baseline) ** 2
            rec_var.append((1 - sa_err / base_err))
            explained_variance_dict[k] = rec_var[i+1] - rec_var[i]
            i += 1

        # Create a DataFrame from the collected data
        explained_variance_df = pd.DataFrame(
            list(explained_variance_dict.items()), 
            columns=['iModulon', 'Explained Variance']
        )
        
        return explained_variance_df