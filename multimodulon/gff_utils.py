"""GFF parsing utilities for MultiModulon."""

from __future__ import annotations

import pandas as pd
from pathlib import Path
from typing import List, Optional, Union, TYPE_CHECKING
import logging

if TYPE_CHECKING:
    from .core import MultiModulon

logger = logging.getLogger(__name__)


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
    if attr_id == 'locus_tag':
        if 'Note' in attributes:
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
        # For gene entries without Note field (W3110 case), use gene name as placeholder
        elif 'gene' in attributes:
            return attributes['gene']
    
    if ignore:
        return None
    else:
        raise ValueError(f"Attribute {attr_id} not found in: {attr_string}")


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
        # Only keep genes that have a valid locus_tag
        DF_gene = DF_gene[DF_gene.locus_tag.notnull()]

        # Sort by start position
        DF_cds = DF_cds.sort_values("start")

        # Extract attribute information
        DF_cds["locus_tag"] = DF_cds.attributes.apply(_get_attr, attr_id="locus_tag", ignore=True)

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
            
            # Store in species data
            species_data._gene_table = df_annot
            
            print(f"  ✓ Created gene table with {len(df_annot)} genes")
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