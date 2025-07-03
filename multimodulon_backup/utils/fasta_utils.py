"""Utilities for working with FASTA files."""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from typing import Dict, Optional


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