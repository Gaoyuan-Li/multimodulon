"""Consolidated utility functions for MultiModulon package."""

from __future__ import annotations

# Standard library imports
import logging
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, TYPE_CHECKING

# Third-party imports
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm.auto import tqdm

# Type checking imports
if TYPE_CHECKING:
    from .core import MultiModulon

# Set up logging
logger = logging.getLogger(__name__)


# ================================================================================
# GFF PARSING UTILITIES
# ================================================================================

def _get_attr(attr_string: str, attr_id: str, ignore: bool = False) -> Optional[str]:
    """
    Extract attribute value from GFF attributes string.
    
    Parameters
    ----------
    attr_string : str
        The attributes string from GFF file
    attr_id : str
        The attribute ID to extract
    ignore : bool
        If True, return None if attribute not found; if False, raise error
        
    Returns
    -------
    str or None
        The attribute value
    """
    attributes = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key] = value
    
    if attr_id in attributes:
        return attributes[attr_id]
    
    # Special handling for locus_tag - try to extract from Note field
    if attr_id == 'locus_tag' and 'Note' in attributes:
        note_value = attributes['Note']
        # Parse "ECK0001:JW4367:b0001" format
        if ':' in note_value:
            parts = note_value.split(':')
            # Look for JW-number first (preferred for W3110)
            for part in parts:
                if part.startswith('JW') and len(part) > 2:
                    return part
            # If no JW-number found, fall back to b-number
            for part in reversed(parts):
                if part.startswith('b') and len(part) > 1:
                    return part
    
    if ignore:
        return None
    else:
        raise ValueError(f"Attribute {attr_id} not found in: {attr_string}")


def _parse_attributes(attr_string: str) -> Dict[str, str]:
    """
    Parse the attributes column of a GFF file.
    
    Parameters
    ----------
    attr_string : str
        The attributes string from the 9th column of a GFF file
        
    Returns
    -------
    Dict[str, str]
        Dictionary of attribute key-value pairs
    """
    attributes = {}
    
    # Split by semicolon
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key] = value
    
    return attributes


def gff2pandas(gff_file: Union[str, List[str]], feature: Union[str, List[str]] = "CDS", 
               index: Optional[str] = None) -> pd.DataFrame:
    """
    Converts GFF file to a Pandas DataFrame with enhanced attribute parsing.
    
    Parameters
    ----------
    gff_file : str or list
        Path(s) to GFF file
    feature: str or list
        Name(s) of features to keep (default = "CDS")
    index : str, optional
        Column or attribute to use as index

    Returns
    -------
    df_gff: pd.DataFrame
        GFF formatted as a DataFrame
    """
    # Argument checking
    if isinstance(gff_file, str):
        gff_file = [gff_file]

    if isinstance(feature, str):
        feature = [feature]

    result = []

    for gff in gff_file:
        with open(gff, "r") as f:
            lines = f.readlines()

        # Get lines to skip
        skiprow = sum([line.startswith("#") for line in lines])

        # Read GFF
        names = [
            "accession",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ]
        DF_gff = pd.read_csv(gff, sep="\t", skiprows=skiprow, names=names, header=None)

        # Filter for CDSs
        DF_cds = DF_gff[DF_gff.feature.isin(feature)]

        # Also filter for genes to get old_locus_tag
        DF_gene = DF_gff[DF_gff.feature == "gene"].reset_index()
        DF_gene["locus_tag"] = DF_gene.attributes.apply(
            _get_attr, attr_id="locus_tag", ignore=True
        )
        DF_gene["old_locus_tag"] = DF_gene.attributes.apply(
            _get_attr, attr_id="old_locus_tag", ignore=True
        )
        DF_gene = DF_gene[["locus_tag", "old_locus_tag"]]
        DF_gene = DF_gene[DF_gene.locus_tag.notnull()]

        # Sort by start position
        DF_cds = DF_cds.sort_values("start")

        # Extract attribute information
        DF_cds["locus_tag"] = DF_cds.attributes.apply(_get_attr, attr_id="locus_tag")

        DF_cds["gene_name"] = DF_cds.attributes.apply(
            _get_attr, attr_id="gene", ignore=True
        )

        DF_cds["gene_product"] = DF_cds.attributes.apply(
            _get_attr, attr_id="product", ignore=True
        )

        DF_cds["ncbi_protein"] = DF_cds.attributes.apply(
            _get_attr, attr_id="protein_id", ignore=True
        )

        # Merge in old_locus_tag
        DF_cds = pd.merge(DF_cds, DF_gene, how="left", on="locus_tag", sort=False)

        result.append(DF_cds)

    DF_gff = pd.concat(result)

    if index:
        if DF_gff[index].duplicated().any():
            logger.debug("Duplicate {} detected. Dropping duplicates.".format(index))
            DF_gff = DF_gff.drop_duplicates(index, keep='first')
        DF_gff.set_index("locus_tag", drop=True, inplace=True)

    return DF_gff


def parse_gff(gff_file: Path) -> pd.DataFrame:
    """
    Parse a GFF file and extract gene information.
    
    Parameters
    ----------
    gff_file : Path
        Path to the GFF file
        
    Returns
    -------
    pd.DataFrame
        DataFrame with gene information indexed by locus_tag
    """
    genes = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            if parts[2] == 'CDS':  # We want protein coding sequences
                # Parse attributes
                attributes = _parse_attributes(parts[8])
                
                # Extract locus_tag
                locus_tag = attributes.get('locus_tag', '')
                if not locus_tag:
                    continue
                
                # Extract gene information
                gene_info = {
                    'locus_tag': locus_tag,
                    'accession': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'gene_name': attributes.get('Name', ''),
                    'old_locus_tag': attributes.get('old_locus_tag', ''),
                    'gene_product': attributes.get('product', '').replace('%2C', ',').replace('%3B', ';'),
                    'ncbi_protein': attributes.get('protein_id', '')
                }
                
                genes.append(gene_info)
    
    # Create DataFrame and set index
    df = pd.DataFrame(genes)
    if not df.empty:
        df = df.set_index('locus_tag')
        # Keep only the requested columns in the specified order
        columns = ['accession', 'start', 'end', 'strand', 'gene_name', 
                   'old_locus_tag', 'gene_product', 'ncbi_protein']
        df = df[columns]
    
    return df


def create_gene_table(multimodulon: 'MultiModulon') -> None:
    """
    Create gene tables for all species by reading GFF files from ref_genome folders.
    
    This method reads GFF files from the ref_genome subfolder in each species directory
    and creates gene tables with columns: accession, start, end, strand, gene_name,
    old_locus_tag, gene_product, ncbi_protein.
    
    The gene tables are stored in each SpeciesData object and can be accessed via
    multimodulon_obj[species_name].gene_table
    
    Parameters
    ----------
    multimodulon : MultiModulon
        The MultiModulon instance
    """
    print("\nCreating gene tables for all species...")
    print("=" * 60)
    
    for species_name, species_data in multimodulon._species_data.items():
        print(f"\nProcessing {species_name}...")
        
        # Find GFF file in ref_genome folder
        ref_genome_path = species_data.data_path / "ref_genome"
        if not ref_genome_path.exists():
            logger.warning(f"No ref_genome folder found for {species_name}, skipping")
            continue
        
        # Look for GFF files
        gff_files = list(ref_genome_path.glob("*.gff"))
        if not gff_files:
            # Also check for .gff3 files
            gff_files = list(ref_genome_path.glob("*.gff3"))
        
        if not gff_files:
            logger.warning(f"No GFF file found in {ref_genome_path}, skipping")
            continue
        
        if len(gff_files) > 1:
            logger.warning(f"Multiple GFF files found for {species_name}, using first: {gff_files[0].name}")
        
        gff_file = gff_files[0]
        logger.info(f"Using GFF file: {gff_file}")
        
        # Parse GFF file using the enhanced parser
        try:
            df_annot = gff2pandas(str(gff_file), index='locus_tag')
            
            # Keep only the specified columns
            keep_cols = ['accession', 'start', 'end', 'strand', 'gene_name', 
                       'old_locus_tag', 'gene_product', 'ncbi_protein']
            
            # Filter for columns that exist
            existing_cols = [col for col in keep_cols if col in df_annot.columns]
            df_annot = df_annot[existing_cols]
            
            # Filter to only keep genes that are in log_tpm matrices
            if species_data.log_tpm is not None:
                # Get genes that are in log_tpm (rows)
                tpm_genes = species_data.log_tpm.index.tolist()
                # Filter gene table to only include these genes
                genes_in_both = [g for g in tpm_genes if g in df_annot.index]
                df_annot = df_annot.loc[genes_in_both]
                
                # Report any genes in log_tpm but not in annotation
                missing_in_annot = set(tpm_genes) - set(df_annot.index)
                if missing_in_annot:
                    logger.warning(f"{species_name}: {len(missing_in_annot)} genes in log_tpm not found in GFF annotation")
                    
                print(f"  ✓ Created gene table with {len(df_annot)} genes (filtered to match log_tpm)")
            else:
                # If no log_tpm data, keep all genes but warn
                logger.warning(f"{species_name}: No log_tpm data found, keeping all {len(df_annot)} genes from GFF")
                print(f"  ✓ Created gene table with {len(df_annot)} genes (no filtering applied)")
            
            # Store in species data
            species_data._gene_table = df_annot
            
            logger.info(f"Gene table for {species_name}: {df_annot.shape}")
            
        except Exception as e:
            logger.error(f"Error processing GFF file for {species_name}: {e}")
            print(f"  ✗ Error: {e}")
    
    print("\n" + "=" * 60)
    print("Gene table creation completed!")
        
    if multimodulon._bbh is not None:
        print(f"\nBBH analysis completed:")
        print(f"  - Total species pairs: {len(multimodulon._bbh) // 2}")
        
        # Print some BBH statistics
        for (sp1, sp2), bbh_df in list(multimodulon._bbh.items())[:3]:
            if not bbh_df.empty:
                print(f"  - {sp1} vs {sp2}: {len(bbh_df)} orthologs")


# ================================================================================
# FASTA UTILITIES
# ================================================================================

def extract_protein_sequences(genome_fasta: Path, gff_file: Path, output_fasta: Path) -> Path:
    """
    Extract protein sequences from genome using GFF annotations.
    
    Parameters
    ----------
    genome_fasta : Path
        Path to genome FASTA file
    gff_file : Path
        Path to GFF annotation file
    output_fasta : Path
        Path for output protein FASTA file
        
    Returns
    -------
    Path
        Path to the created protein FASTA file
    """
    # Read genome sequence
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    
    # Parse GFF to get CDS features
    protein_records = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            if parts[2] == 'CDS':
                # Parse attributes
                attributes = {}
                for attr in parts[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value
                
                # Get required information
                accession = parts[0]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                strand = parts[6]
                locus_tag = attributes.get('locus_tag', '')
                protein_id = attributes.get('protein_id', '')
                product = attributes.get('product', '').replace('%2C', ',').replace('%3B', ';')
                
                if not locus_tag or accession not in genome_dict:
                    continue
                
                # Extract sequence
                seq = genome_dict[accession].seq[start:end]
                
                if strand == '-':
                    seq = seq.reverse_complement()
                
                # Translate to protein
                try:
                    # Trim sequence to multiple of 3 to avoid partial codon warning
                    seq_len = len(seq)
                    if seq_len % 3 != 0:
                        seq = seq[:-(seq_len % 3)]
                    
                    protein_seq = seq.translate(to_stop=True)
                    
                    # Create SeqRecord
                    # Build description carefully to avoid empty strings
                    desc_parts = []
                    if protein_id:
                        desc_parts.append(protein_id)
                    if product:
                        desc_parts.append(product)
                    
                    record = SeqRecord(
                        protein_seq,
                        id=locus_tag,
                        description=" ".join(desc_parts) if desc_parts else locus_tag
                    )
                    protein_records.append(record)
                    
                except Exception as e:
                    # Skip sequences that can't be translated
                    continue
    
    # Write protein sequences
    if protein_records:
        SeqIO.write(protein_records, output_fasta, "fasta")
    
    return output_fasta


# ================================================================================
# BIDIRECTIONAL BEST HITS (BBH) UTILITIES
# ================================================================================

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