#!/usr/bin/env python3
"""
Example usage of the MultiModulon package.
"""

from multimodulon import MultiModulon
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def main():
    # Path to the data folder
    data_folder = '/home/gaoyuan/multimodulon_development/imminer_2_pseudomonas'
    
    print("Initializing MultiModulon...")
    print("=" * 50)
    
    # Initialize MultiModulon object
    # Note: Set run_bbh=False to skip BBH analysis for faster loading
    multi_modulon = MultiModulon(data_folder, run_bbh=False)
    
    # Print summary
    print("\n")
    multi_modulon.summary()
    
    # Access specific species data
    print("\n" + "=" * 50)
    print("Accessing Pseudomonas_fluorescens data:")
    print("=" * 50)
    
    pf_data = multi_modulon['Pseudomonas_fluorescens']
    
    # Access expression data
    print(f"\nLog TPM matrix shape: {pf_data.log_tpm.shape}")
    print(f"First 5 genes: {list(pf_data.log_tpm.index[:5])}")
    print(f"First 5 samples: {list(pf_data.log_tpm.columns[:5])}")
    
    # Access normalized expression data
    print(f"\nNormalized expression matrix (X) shape: {pf_data.X.shape}")
    
    # Access sample metadata
    print(f"\nSample sheet shape: {pf_data.sample_sheet.shape}")
    print(f"Sample sheet columns: {list(pf_data.sample_sheet.columns[:10])}")
    
    # Access gene annotations
    print(f"\nGene table shape: {pf_data.gene_table.shape}")
    print(f"Gene table columns: {list(pf_data.gene_table.columns)}")
    print("\nFirst 5 genes from gene table:")
    print(pf_data.gene_table.head())
    
    # Run BBH analysis (if needed)
    print("\n" + "=" * 50)
    print("Running BBH analysis (this may take a while)...")
    print("=" * 50)
    
    # Uncomment to run BBH analysis
    # multi_modulon._run_bbh_analysis()
    
    # Save BBH results
    # multi_modulon.save_bbh('bbh_results')
    
    # Get orthologs between species (after BBH is run)
    # orthologs = multi_modulon.get_orthologs('Pseudomonas_fluorescens', 'Pseudomonas_putida')
    # print(f"\nFound {len(orthologs)} orthologs between P. fluorescens and P. putida")
    
    print("\nExample completed!")


if __name__ == "__main__":
    main()