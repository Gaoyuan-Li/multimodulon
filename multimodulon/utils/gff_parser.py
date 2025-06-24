"""GFF file parser for extracting gene information."""

import pandas as pd
from typing import Dict, List, Optional
import re
from pathlib import Path


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