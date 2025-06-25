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
    multi_modulon = MultiModulon(input_data_path)
    
    # Print summary
    print("\n")
    multi_modulon.summary()
    
    # Access specific strain data
    print("\n" + "=" * 60)
    print("Accessing strain data:")
    print("=" * 60)
    
    # Get available strains
    available_strains = multi_modulon.species
    print(f"Available strains: {available_strains}")
    
    if available_strains:
        # Access first strain data
        first_strain = available_strains[0]
        strain_data = multi_modulon[first_strain]
        
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
    
    # This will create BBH files for all strain pairs
    multi_modulon.generate_BBH("../imminer_2_industrial_strain/Output_BBH")
    
    # Align genes across all strains
    print("\n" + "=" * 60)
    print("Aligning genes across strains...")
    print("=" * 60)
    
    # This will create the combined gene database and aligned expression matrices
    combined_gene_db = multi_modulon.align_genes("../imminer_2_industrial_strain/Output_Gene_Info")
    
    print(f"\nCombined gene database shape: {combined_gene_db.shape}")
    print(f"Columns: {list(combined_gene_db.columns)}")
    
    # Access aligned expression matrices
    print("\n" + "=" * 60)
    print("Accessing aligned expression matrices:")
    print("=" * 60)
    
    for strain_name in available_strains:
        strain_data = multi_modulon[strain_name]
        print(f"{strain_name} aligned X matrix shape: {strain_data.X.shape}")
        print(f"  Index (first 5): {list(strain_data.X.index[:5])}")
    
    print("\nExample completed successfully!")
    print("You can now access:")
    print("- multi_modulon['strain_name'].X for aligned expression matrices")
    print("- Combined gene database saved to Output_Gene_Info/combined_gene_db.csv")
    print("- BBH files saved to Output_BBH/")


if __name__ == "__main__":
    main()