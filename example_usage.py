#!/usr/bin/env python3
"""
Example usage of the updated MultiModulon package with CSV support.
"""

from multimodulon import MultiModulon
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def main():
    # Path to the Input_Data folder
    input_data_path = '../imminer_2_industrial_strain/Input_Data'
    
    print("Initializing MultiModulon with new CSV-based structure...")
    print("=" * 60)
    
    # Initialize MultiModulon object with the Input_Data path
    multiModulon = MultiModulon(input_data_path)
    
    # Print summary
    print("\n")
    multiModulon.summary()
    
    # Access specific strain data
    print("\n" + "=" * 60)
    print("Accessing strain data:")
    print("=" * 60)
    
    # Get available strains
    available_strains = multiModulon.species
    print(f"Available strains: {available_strains}")
    
    if available_strains:
        # Access first strain data
        first_strain = available_strains[0]
        strain_data = multiModulon[first_strain]
        
        print(f"\n{first_strain} data:")
        print(f"  - Log TPM matrix shape: {strain_data.log_tpm.shape}")
        print(f"  - Log TPM normalized matrix shape: {strain_data.log_tpm_norm.shape}")
        print(f"  - X matrix shape (alias for log_tpm_norm): {strain_data.X.shape}")
        print(f"  - Sample sheet shape: {strain_data.sample_sheet.shape}")
        print(f"  - First 5 genes: {list(strain_data.log_tpm.index[:5])}")
        print(f"  - First 5 samples: {list(strain_data.log_tpm.columns[:5])}")
    
    # Generate BBH files
    print("\n" + "=" * 60)
    print("Generating BBH files...")
    print("=" * 60)
    
    # This will create BBH files for all strain pairs using existing protein.faa files
    # The BBH files will use locus_tags (e.g., b0001, EcolC_0001) in the gene and subject columns
    # Use threads parameter to speed up BLAST computation (default is 1)
    multiModulon.generate_BBH("../imminer_2_industrial_strain/Output_BBH_test", threads=8)
    
    # Align genes across all strains
    print("\n" + "=" * 60)
    print("Aligning genes across strains...")
    print("=" * 60)
    
    # This will create the combined gene database and aligned expression matrices
    combined_gene_db = multiModulon.align_genes(
        input_bbh_dir="../imminer_2_industrial_strain/Output_BBH_test",
        output_dir="../imminer_2_industrial_strain/Output_Gene_Info"
    )
    
    print(f"\nCombined gene database shape: {combined_gene_db.shape}")
    print(f"Columns: {list(combined_gene_db.columns)}")
    
    # Create gene tables from GFF files
    print("\n" + "=" * 60)
    print("Creating gene tables from GFF files...")
    print("=" * 60)
    
    # This will read GFF files from ref_genome folder in each species directory
    multiModulon.create_gene_table()
    
    # Access gene tables
    print("\n" + "=" * 60)
    print("Accessing gene tables:")
    print("=" * 60)
    
    for strain_name in available_strains:
        strain_data = multiModulon[strain_name]
        if hasattr(strain_data, '_gene_table') and strain_data._gene_table is not None:
            print(f"\n{strain_name} gene table:")
            print(f"  - Shape: {strain_data.gene_table.shape}")
            print(f"  - Columns: {list(strain_data.gene_table.columns)}")
            print(f"  - First 3 genes:")
            print(strain_data.gene_table.head(3))
    
    # Access aligned expression matrices
    print("\n" + "=" * 60)
    print("Accessing aligned expression matrices:")
    print("=" * 60)
    
    for strain_name in available_strains:
        strain_data = multiModulon[strain_name]
        print(f"{strain_name} aligned X matrix shape: {strain_data.X.shape}")
        print(f"  Index (first 5): {list(strain_data.X.index[:5])}")
    
    print("\nExample completed successfully!")
    print("You can now access:")
    print("- multiModulon['strain_name'].gene_table for gene annotations")
    print("- multiModulon['strain_name'].X for aligned expression matrices")
    print("- Combined gene database saved to Output_Gene_Info/combined_gene_db.csv")
    print("- BBH files saved to Output_BBH/")


if __name__ == "__main__":
    main()