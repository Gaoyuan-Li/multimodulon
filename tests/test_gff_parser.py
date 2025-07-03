"""Tests for GFF parser functionality."""

import pytest
import tempfile
from pathlib import Path
import pandas as pd

from multimodulon.utils import parse_gff, _parse_attributes


class TestGFFParser:
    """Test GFF parsing functionality."""
    
    def test_parse_attributes(self):
        """Test parsing of GFF attributes."""
        attr_string = "ID=cds-WP_096237022.1;Parent=gene-CPH89_RS00005;Name=WP_096237022.1;locus_tag=CPH89_RS00005"
        attrs = _parse_attributes(attr_string)
        
        assert attrs['ID'] == 'cds-WP_096237022.1'
        assert attrs['Parent'] == 'gene-CPH89_RS00005'
        assert attrs['Name'] == 'WP_096237022.1'
        assert attrs['locus_tag'] == 'CPH89_RS00005'
    
    def test_parse_gff_basic(self):
        """Test basic GFF parsing."""
        # Create a temporary GFF file
        gff_content = """##gff-version 3
NZ_LT907842.1	RefSeq	gene	1	514	.	-	.	ID=gene-CPH89_RS00005;Name=CPH89_RS00005;locus_tag=CPH89_RS00005
NZ_LT907842.1	RefSeq	CDS	1	514	.	-	0	ID=cds-WP_096237022.1;Parent=gene-CPH89_RS00005;Name=WP_096237022.1;locus_tag=CPH89_RS00005;product=transporter;protein_id=WP_096237022.1;old_locus_tag=OLD_001
NZ_LT907842.1	RefSeq	gene	756	1181	.	-	.	ID=gene-CPH89_RS00010;Name=CPH89_RS00010;locus_tag=CPH89_RS00010
NZ_LT907842.1	RefSeq	CDS	756	1181	.	-	0	ID=cds-WP_053257077.1;Parent=gene-CPH89_RS00010;Name=WP_053257077.1;locus_tag=CPH89_RS00010;product=acetyltransferase;protein_id=WP_053257077.1
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write(gff_content)
            temp_path = Path(f.name)
        
        try:
            # Parse the GFF file
            gene_table = parse_gff(temp_path)
            
            # Check the results
            assert len(gene_table) == 2
            assert 'CPH89_RS00005' in gene_table.index
            assert 'CPH89_RS00010' in gene_table.index
            
            # Check columns
            expected_columns = ['accession', 'start', 'end', 'strand', 'gene_name', 
                              'old_locus_tag', 'gene_product', 'ncbi_protein']
            assert list(gene_table.columns) == expected_columns
            
            # Check specific values
            gene1 = gene_table.loc['CPH89_RS00005']
            assert gene1['accession'] == 'NZ_LT907842.1'
            assert gene1['start'] == 1
            assert gene1['end'] == 514
            assert gene1['strand'] == '-'
            assert gene1['gene_name'] == 'WP_096237022.1'
            assert gene1['old_locus_tag'] == 'OLD_001'
            assert gene1['gene_product'] == 'transporter'
            assert gene1['ncbi_protein'] == 'WP_096237022.1'
            
        finally:
            # Clean up
            temp_path.unlink()
    
    def test_parse_gff_empty(self):
        """Test parsing empty GFF file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("##gff-version 3\n")
            temp_path = Path(f.name)
        
        try:
            gene_table = parse_gff(temp_path)
            assert len(gene_table) == 0
            assert isinstance(gene_table, pd.DataFrame)
        finally:
            temp_path.unlink()
    
    def test_parse_gff_special_characters(self):
        """Test parsing GFF with special characters in product names."""
        gff_content = """##gff-version 3
NZ_LT907842.1	RefSeq	CDS	1	514	.	-	0	ID=cds-WP_096237022.1;locus_tag=CPH89_RS00005;product=ABC transporter%2C ATP-binding protein%3B membrane protein;protein_id=WP_096237022.1
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write(gff_content)
            temp_path = Path(f.name)
        
        try:
            gene_table = parse_gff(temp_path)
            
            # Check that special characters are properly decoded
            assert len(gene_table) == 1
            gene = gene_table.iloc[0]
            assert gene['gene_product'] == 'ABC transporter, ATP-binding protein; membrane protein'
            
        finally:
            temp_path.unlink()