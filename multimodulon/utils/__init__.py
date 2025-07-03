"""Utility functions for multimodulon package."""

# Import from the consolidated utils module
from ..utils import (
    parse_gff,
    gff2pandas,
    create_gene_table,
    BBHAnalyzer,
    extract_protein_sequences,
    _get_attr,
    _parse_attributes
)

__all__ = [
    'parse_gff',
    'gff2pandas', 
    'create_gene_table',
    'BBHAnalyzer',
    'extract_protein_sequences',
    '_get_attr',
    '_parse_attributes'
]