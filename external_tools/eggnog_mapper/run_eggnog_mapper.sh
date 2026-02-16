#!/bin/bash

# run_eggnog_mapper.sh
# Script to run eggNOG-mapper on multiple species from imminer_2_industrial_strain dataset
# Usage: ./run_eggnog_mapper.sh <input_data_path> <output_dir_path> <tax_scope>
# Example: ./run_eggnog_mapper.sh ../imminer_2_industrial_strain/Input_Data ../imminer_2_industrial_strain/Output_eggnog_mapper Gammaproteobacteria

# Check if correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_data_path> <output_dir_path> <tax_scope>"
    echo "Example: $0 ../imminer_2_industrial_strain/Input_Data ../imminer_2_industrial_strain/Output_eggnog_mapper Gammaproteobacteria"
    exit 1
fi

# Set input and output paths
INPUT_DATA_PATH="$1"
OUTPUT_DIR="$2"
TAX_SCOPE="$3"

# Check if input path exists
if [ ! -d "$INPUT_DATA_PATH" ]; then
    echo "Error: Input data path does not exist: $INPUT_DATA_PATH"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Get species list from subfolder names in Input_Data
SPECIES=()
for dir in "$INPUT_DATA_PATH"/*; do
    if [ -d "$dir" ]; then
        SPECIES+=("$(basename "$dir")")
    fi
done

echo "Found species: ${SPECIES[@]}"

# Process each species
for species in "${SPECIES[@]}"; do
    echo "Processing species: $species"
    
    # Set input faa file path
    FAA_FILE="$INPUT_DATA_PATH/$species/ref_genome/protein.faa"
    
    # Check if faa file exists
    if [ ! -f "$FAA_FILE" ]; then
        echo "Warning: Protein file not found for $species: $FAA_FILE"
        echo "Skipping $species..."
        continue
    fi
    
    # Create output directory for this species
    SPECIES_OUTPUT_DIR="$OUTPUT_DIR/$species"
    mkdir -p "$SPECIES_OUTPUT_DIR"
    
    # Run eggNOG-mapper
    echo "Running eggNOG-mapper for $species..."
    python /mnt/bop/eggnog-mapper/emapper.py \
        -o "$species" \
        --tax_scope "$TAX_SCOPE" \
        --tax_scope_mode Bacteria \
        -i "$FAA_FILE" \
        --output_dir "$SPECIES_OUTPUT_DIR" \
        --cpu 16
    
    # Check if command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed $species"
    else
        echo "Error processing $species"
    fi
    
    echo "----------------------------------------"
done

echo "All species processed!"