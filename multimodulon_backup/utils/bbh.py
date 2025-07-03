"""Bidirectional Best Hits (BBH) implementation for ortholog detection."""

import subprocess
import tempfile
from pathlib import Path
import pandas as pd
from typing import Dict, Tuple, Optional
import os
from tqdm.auto import tqdm
import logging

logger = logging.getLogger(__name__)


class BBHAnalyzer:
    """
    Performs Bidirectional Best Hits analysis between species.
    """
    
    def __init__(self, docker_image: str = "quay.io/biocontainers/blast:2.16.0--h66d330f_5"):
        """
        Initialize BBH analyzer.
        
        Parameters
        ----------
        docker_image : str
            Docker image to use for BLAST
        """
        self.docker_image = docker_image
        self._check_docker()
    
    def _check_docker(self):
        """Check if Docker is available."""
        try:
            subprocess.run(["docker", "--version"], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("Docker not found. BBH analysis may fail.")
    
    def run_bbh_all_pairs(self, species_data: Dict[str, Dict]) -> Dict[Tuple[str, str], pd.DataFrame]:
        """
        Run BBH analysis for all pairs of species.
        
        Parameters
        ----------
        species_data : Dict[str, Dict]
            Dictionary with species names as keys and data dictionaries containing
            'fasta_path' with path to protein FASTA file
            
        Returns
        -------
        Dict[Tuple[str, str], pd.DataFrame]
            Dictionary with species pairs as keys and BBH results as values
        """
        bbh_results = {}
        species_list = list(species_data.keys())
        
        # Create all pairs
        pairs = []
        for i in range(len(species_list)):
            for j in range(i + 1, len(species_list)):
                pairs.append((species_list[i], species_list[j]))
        
        # Run BBH for each pair
        for sp1, sp2 in tqdm(pairs, desc="Running BBH analysis"):
            logger.info(f"Running BBH between {sp1} and {sp2}")
            
            try:
                bbh_df = self._run_bbh_pair(
                    species_data[sp1]['fasta_path'],
                    species_data[sp2]['fasta_path'],
                    sp1, sp2
                )
                
                # Store results bidirectionally
                bbh_results[(sp1, sp2)] = bbh_df
                bbh_results[(sp2, sp1)] = bbh_df[['query_id', 'subject_id', 'evalue', 'bitscore']]
                bbh_results[(sp2, sp1)].columns = ['subject_id', 'query_id', 'evalue', 'bitscore']
                
            except Exception as e:
                logger.error(f"Failed to run BBH between {sp1} and {sp2}: {str(e)}")
                bbh_results[(sp1, sp2)] = pd.DataFrame()
                bbh_results[(sp2, sp1)] = pd.DataFrame()
        
        return bbh_results
    
    def _run_bbh_pair(self, fasta1: Path, fasta2: Path, sp1: str, sp2: str) -> pd.DataFrame:
        """
        Run BBH analysis between two species.
        
        Parameters
        ----------
        fasta1, fasta2 : Path
            Paths to protein FASTA files
        sp1, sp2 : str
            Species names
            
        Returns
        -------
        pd.DataFrame
            BBH results with columns: query_id, subject_id, evalue, bitscore
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Copy FASTA files to temp directory
            tmp_fasta1 = tmpdir / f"{sp1}.faa"
            tmp_fasta2 = tmpdir / f"{sp2}.faa"
            
            subprocess.run(f"cp {fasta1} {tmp_fasta1}", shell=True, check=True)
            subprocess.run(f"cp {fasta2} {tmp_fasta2}", shell=True, check=True)
            
            # Run BLAST in both directions
            blast1_out = tmpdir / f"{sp1}_vs_{sp2}.blast"
            blast2_out = tmpdir / f"{sp2}_vs_{sp1}.blast"
            
            # First direction: sp1 query vs sp2 database
            self._run_blast(tmp_fasta1, tmp_fasta2, blast1_out, tmpdir)
            
            # Second direction: sp2 query vs sp1 database
            self._run_blast(tmp_fasta2, tmp_fasta1, blast2_out, tmpdir)
            
            # Parse BLAST results
            hits1 = self._parse_blast_output(blast1_out)
            hits2 = self._parse_blast_output(blast2_out)
            
            # Find bidirectional best hits
            bbh = self._find_bbh(hits1, hits2)
            
            return bbh
    
    def _run_blast(self, query: Path, subject: Path, output: Path, tmpdir: Path):
        """Run BLAST using Docker."""
        # Create BLAST database
        db_name = subject.stem
        
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{tmpdir}:/data",
            self.docker_image,
            "makeblastdb",
            "-in", f"/data/{subject.name}",
            "-dbtype", "prot",
            "-out", f"/data/{db_name}"
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        
        # Run BLAST
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{tmpdir}:/data",
            self.docker_image,
            "blastp",
            "-query", f"/data/{query.name}",
            "-db", f"/data/{db_name}",
            "-out", f"/data/{output.name}",
            "-outfmt", "6 qseqid sseqid evalue bitscore",
            "-max_target_seqs", "1",
            "-evalue", "1e-5"
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
    
    def _parse_blast_output(self, blast_file: Path) -> pd.DataFrame:
        """Parse BLAST output file."""
        if not blast_file.exists() or blast_file.stat().st_size == 0:
            return pd.DataFrame(columns=['query_id', 'subject_id', 'evalue', 'bitscore'])
        
        df = pd.read_csv(
            blast_file,
            sep='\t',
            names=['query_id', 'subject_id', 'evalue', 'bitscore'],
            dtype={'query_id': str, 'subject_id': str, 'evalue': float, 'bitscore': float}
        )
        
        # Keep only best hit per query
        df = df.sort_values(['query_id', 'bitscore'], ascending=[True, False])
        df = df.groupby('query_id').first().reset_index()
        
        return df
    
    def _find_bbh(self, hits1: pd.DataFrame, hits2: pd.DataFrame) -> pd.DataFrame:
        """Find bidirectional best hits."""
        if hits1.empty or hits2.empty:
            return pd.DataFrame(columns=['query_id', 'subject_id', 'evalue', 'bitscore'])
        
        # Create dictionaries for fast lookup
        forward_hits = dict(zip(hits1['query_id'], hits1['subject_id']))
        reverse_hits = dict(zip(hits2['query_id'], hits2['subject_id']))
        
        # Find BBH
        bbh_pairs = []
        for query, subject in forward_hits.items():
            if subject in reverse_hits and reverse_hits[subject] == query:
                # This is a BBH pair
                row = hits1[hits1['query_id'] == query].iloc[0]
                bbh_pairs.append(row)
        
        if bbh_pairs:
            bbh_df = pd.DataFrame(bbh_pairs)
        else:
            bbh_df = pd.DataFrame(columns=['query_id', 'subject_id', 'evalue', 'bitscore'])
        
        return bbh_df