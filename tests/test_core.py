"""Tests for core MultiModulon functionality."""

import pytest
import tempfile
from pathlib import Path
import pandas as pd
import numpy as np

from multimodulon.core import SpeciesData, MultiModulon


class TestSpeciesData:
    """Test SpeciesData functionality."""
    
    def create_test_data(self, tmpdir: Path) -> Path:
        """Create test data structure."""
        species_dir = tmpdir / "Test_species"
        species_dir.mkdir()
        
        # Create expression matrices directory
        expr_dir = species_dir / "expression_matrices"
        expr_dir.mkdir()
        
        # Create test log_tpm.tsv
        log_tpm_data = pd.DataFrame({
            'Sample1': [1.0, 2.0, 3.0],
            'Sample2': [1.5, 2.5, 3.5],
            'Sample3': [2.0, 3.0, 4.0]
        }, index=['gene-CPH89_RS00005', 'gene-CPH89_RS00010', 'gene-CPH89_RS00015'])
        log_tpm_data.to_csv(expr_dir / 'log_tpm.tsv', sep='\t')
        
        # Create test log_tpm_norm.csv
        log_tpm_norm_data = log_tpm_data - log_tpm_data.mean()
        log_tpm_norm_data.to_csv(expr_dir / 'log_tpm_norm.csv')
        
        # Create samplesheet directory
        sample_dir = species_dir / "samplesheet"
        sample_dir.mkdir()
        
        # Create test samplesheet.csv
        sample_data = pd.DataFrame({
            'sample': ['Sample1', 'Sample2', 'Sample3'],
            'condition': ['Control', 'Treatment', 'Treatment'],
            'replicate': [1, 1, 2]
        })
        sample_data.to_csv(sample_dir / 'samplesheet.csv', index=False)
        
        # Create ref_genome directory
        genome_dir = species_dir / "ref_genome"
        genome_dir.mkdir()
        
        # Create test GFF file
        gff_content = """##gff-version 3
NZ_TEST.1	RefSeq	CDS	1	514	.	-	0	ID=cds-WP_001;locus_tag=CPH89_RS00005;product=test protein 1;protein_id=WP_001
NZ_TEST.1	RefSeq	CDS	756	1181	.	-	0	ID=cds-WP_002;locus_tag=CPH89_RS00010;product=test protein 2;protein_id=WP_002
NZ_TEST.1	RefSeq	CDS	1500	2000	.	+	0	ID=cds-WP_003;locus_tag=CPH89_RS00015;product=test protein 3;protein_id=WP_003
"""
        with open(genome_dir / 'genomic.gff', 'w') as f:
            f.write(gff_content)
        
        return species_dir
    
    def test_species_data_loading(self):
        """Test loading species data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            species_dir = self.create_test_data(tmpdir)
            
            # Create SpeciesData object
            species_data = SpeciesData('Test_species', species_dir)
            
            # Test log_tpm loading
            assert species_data.log_tpm.shape == (3, 3)
            assert list(species_data.log_tpm.columns) == ['Sample1', 'Sample2', 'Sample3']
            assert list(species_data.log_tpm.index) == ['gene-CPH89_RS00005', 'gene-CPH89_RS00010', 'gene-CPH89_RS00015']
            
            # Test X (normalized) loading
            assert species_data.X.shape == (3, 3)
            assert np.allclose(species_data.X.mean(), 0, atol=1e-10)
            
            # Test sample_sheet loading
            assert species_data.sample_sheet.shape == (3, 2)  # 3 samples, 2 columns (condition, replicate)
            assert 'Sample1' in species_data.sample_sheet.index
            
            # Test gene_table loading
            assert species_data.gene_table.shape == (3, 8)
            assert 'CPH89_RS00005' in species_data.gene_table.index
            assert species_data.gene_table.loc['CPH89_RS00005', 'gene_product'] == 'test protein 1'
    
    def test_data_validation(self):
        """Test data validation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            species_dir = self.create_test_data(tmpdir)
            
            species_data = SpeciesData('Test_species', species_dir)
            
            # Should pass validation
            assert species_data.validate_data() == True


class TestMultiModulon:
    """Test MultiModulon functionality."""
    
    def create_multi_species_data(self, tmpdir: Path) -> Path:
        """Create test data for multiple species."""
        # Create first species
        sp1_dir = tmpdir / "Species_A"
        sp1_dir.mkdir()
        
        # Expression matrices
        expr_dir1 = sp1_dir / "expression_matrices"
        expr_dir1.mkdir()
        
        log_tpm1 = pd.DataFrame({
            'Sample1': [1.0, 2.0],
            'Sample2': [1.5, 2.5]
        }, index=['gene-GeneA1', 'gene-GeneA2'])
        log_tpm1.to_csv(expr_dir1 / 'log_tpm.tsv', sep='\t')
        log_tpm1.to_csv(expr_dir1 / 'log_tpm_norm.csv')
        
        # Samplesheet
        sample_dir1 = sp1_dir / "samplesheet"
        sample_dir1.mkdir()
        pd.DataFrame({'sample': ['Sample1', 'Sample2']}).to_csv(sample_dir1 / 'samplesheet.csv', index=False)
        
        # Ref genome
        genome_dir1 = sp1_dir / "ref_genome"
        genome_dir1.mkdir()
        
        gff1 = """##gff-version 3
Chr1	RefSeq	CDS	1	100	.	+	0	ID=cds-1;locus_tag=GeneA1;product=protein A1;protein_id=PA1
Chr1	RefSeq	CDS	200	300	.	+	0	ID=cds-2;locus_tag=GeneA2;product=protein A2;protein_id=PA2
"""
        with open(genome_dir1 / 'genomic.gff', 'w') as f:
            f.write(gff1)
        
        # Create second species
        sp2_dir = tmpdir / "Species_B"
        sp2_dir.mkdir()
        
        # Expression matrices
        expr_dir2 = sp2_dir / "expression_matrices"
        expr_dir2.mkdir()
        
        log_tpm2 = pd.DataFrame({
            'Sample3': [3.0, 4.0],
            'Sample4': [3.5, 4.5]
        }, index=['gene-GeneB1', 'gene-GeneB2'])
        log_tpm2.to_csv(expr_dir2 / 'log_tpm.tsv', sep='\t')
        log_tpm2.to_csv(expr_dir2 / 'log_tpm_norm.csv')
        
        # Samplesheet
        sample_dir2 = sp2_dir / "samplesheet"
        sample_dir2.mkdir()
        pd.DataFrame({'sample': ['Sample3', 'Sample4']}).to_csv(sample_dir2 / 'samplesheet.csv', index=False)
        
        # Ref genome
        genome_dir2 = sp2_dir / "ref_genome"
        genome_dir2.mkdir()
        
        gff2 = """##gff-version 3
Chr1	RefSeq	CDS	1	100	.	+	0	ID=cds-1;locus_tag=GeneB1;product=protein B1;protein_id=PB1
Chr1	RefSeq	CDS	200	300	.	+	0	ID=cds-2;locus_tag=GeneB2;product=protein B2;protein_id=PB2
"""
        with open(genome_dir2 / 'genomic.gff', 'w') as f:
            f.write(gff2)
        
        return tmpdir
    
    def test_multiModulon_initialization(self):
        """Test MultiModulon initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            data_dir = self.create_multi_species_data(tmpdir)
            
            # Initialize without BBH
            multi = MultiModulon(str(data_dir), run_bbh=False)
            
            # Check species loaded
            assert len(multi.species) == 2
            assert 'Species_A' in multi.species
            assert 'Species_B' in multi.species
            
            # Access species data
            sp_a = multi['Species_A']
            assert sp_a.log_tpm.shape == (2, 2)
            
            sp_b = multi['Species_B']
            assert sp_b.log_tpm.shape == (2, 2)
    
    def test_multiModulon_invalid_species(self):
        """Test accessing invalid species."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            data_dir = self.create_multi_species_data(tmpdir)
            
            multi = MultiModulon(str(data_dir), run_bbh=False)
            
            with pytest.raises(KeyError) as exc_info:
                _ = multi['Invalid_Species']
            
            assert 'Invalid_Species' in str(exc_info.value)
    
    def test_summary_output(self, capsys):
        """Test summary output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            data_dir = self.create_multi_species_data(tmpdir)
            
            multi = MultiModulon(str(data_dir), run_bbh=False)
            multi.summary()
            
            captured = capsys.readouterr()
            assert 'MultiModulon Summary' in captured.out
            assert 'Species_A' in captured.out
            assert 'Species_B' in captured.out
            assert 'Genes: 2' in captured.out
            assert 'Samples: 2' in captured.out