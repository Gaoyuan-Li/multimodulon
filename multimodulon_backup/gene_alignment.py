"""Gene alignment and BBH analysis functions for MultiModulon."""

from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .core import MultiModulon

logger = logging.getLogger(__name__)


class UnionFind:
    """Union-Find data structure for gene grouping."""
    
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


def generate_BBH(multimodulon: 'MultiModulon', output_path: str = "Output_BBH", threads: int = 1):
    """
    Generate BBH files using existing protein.faa files from each strain.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        The MultiModulon instance
    output_path : str
        Path to save BBH results (default: "Output_BBH")
    threads : int
        Number of threads to use for BLAST (default: 1)
    """
    output_dir = Path(output_path)
    output_dir.mkdir(exist_ok=True)
    
    species_list = list(multimodulon._species_data.keys())
    logger.info(f"Generating BBH for {len(species_list)} species: {species_list}")
    
    # Use existing protein.faa files and create mapping from protein_id to locus_tag
    species_fasta_paths = {}
    protein_to_locus_maps = {}
    
    for species_name in species_list:
        species_data = multimodulon._species_data[species_name]
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
        protein_to_locus = _create_protein_to_locus_mapping(multimodulon, gff_file, species_name)
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
                blast_df = _run_blast_comparison_with_locus(
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


def align_genes(multimodulon: 'MultiModulon', input_bbh_dir: str = "Output_BBH", 
                output_dir: str = "Output_Gene_Info", 
                reference_order: Optional[List[str]] = None, 
                bbh_threshold: Optional[float] = None) -> pd.DataFrame:
    """
    Align genes across all species using Union-Find algorithm to create combined gene database.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        The MultiModulon instance
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
        
    Returns
    -------
    pd.DataFrame
        Combined gene database
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    bbh_dir = Path(input_bbh_dir)
    if not bbh_dir.exists():
        raise FileNotFoundError(f"BBH directory not found: {input_bbh_dir}")
    
    # Get list of strains
    if reference_order:
        # Validate that all strains in reference_order exist in species_data
        available_strains = set(multimodulon._species_data.keys())
        for strain in reference_order:
            if strain not in available_strains:
                raise ValueError(f"Strain '{strain}' in reference_order not found in loaded species data. Available strains: {list(available_strains)}")
        strains = reference_order
    else:
        strains = list(multimodulon._species_data.keys())
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
            species_data = multimodulon._species_data[strain]
            genes = species_data.log_tpm.index.tolist()
            strain_genes[strain].update(genes)
            logger.info(f"Using {len(genes)} genes from expression data for {strain}")
    
    # Step 2: Initialize Union-Find structure
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
    multimodulon._combined_gene_db = df
    
    logger.info(f"Gene alignment completed. Combined gene database saved to {output_file}")
    logger.info(f"Total gene groups: {len(df)}")
    
    # Create aligned expression matrices
    _create_aligned_expression_matrices(multimodulon, df, strains)
    
    return df


def _create_protein_to_locus_mapping(multimodulon: 'MultiModulon', gff_file: Path, species_name: str) -> dict:
    """
    Create mapping from protein_id (NCBI_GP) to locus_tag from GFF file.
    Handles different GFF formats including W3110 which uses Note field instead of locus_tag.
    
    For W3110 strain, extracts JW-numbers (e.g., JW4367) from Note field format "ECK0001:JW4367:b0001".
    Issues warnings when fallback extraction methods are used.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        The MultiModulon instance
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
                    if warning_key not in getattr(multimodulon, '_warned_extractions', set()):
                        print(f"WARNING: '{species_name}' - No locus_tag found in GFF file, extracting from Note field")
                        if not hasattr(multimodulon, '_warned_extractions'):
                            multimodulon._warned_extractions = set()
                        multimodulon._warned_extractions.add(warning_key)
                    
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
                    if warning_key not in getattr(multimodulon, '_warned_extractions', set()):
                        print(f"WARNING: Species '{species_name}' - No locus_tag or Note field found in GFF file, using gene name as fallback")
                        if not hasattr(multimodulon, '_warned_extractions'):
                            multimodulon._warned_extractions = set()
                        multimodulon._warned_extractions.add(warning_key)
                    locus_tag = attr_dict['gene']
                
                # If still no locus_tag, use protein_id itself as last resort
                if not locus_tag and protein_id:
                    warning_key = f"{species_name}_protein_id_fallback"
                    if warning_key not in getattr(multimodulon, '_warned_extractions', set()):
                        print(f"WARNING: Species '{species_name}' - No identifiable locus_tag in GFF file, using protein_id as last resort")
                        if not hasattr(multimodulon, '_warned_extractions'):
                            multimodulon._warned_extractions = set()
                        multimodulon._warned_extractions.add(warning_key)
                    locus_tag = protein_id
                
                if protein_id and locus_tag:
                    protein_to_locus[protein_id] = locus_tag
    
    return protein_to_locus


def _run_blast_comparison_with_locus(fasta1: Path, fasta2: Path, sp1: str, sp2: str,
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


def _create_aligned_expression_matrices(multimodulon: 'MultiModulon', combined_gene_db: pd.DataFrame, strains: List[str]):
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
        species_data = multimodulon._species_data[strain]
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
            multimodulon._species_data[strain]._X = new_df
            logger.info(f"Created aligned expression matrix for {strain}: {new_df.shape}")
        else:
            logger.warning(f"No expression data found for {strain}")