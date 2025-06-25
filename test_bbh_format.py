#!/usr/bin/env python3
"""
Test script to verify BBH output format matches the expected format.
"""

import pandas as pd
import sys

def compare_bbh_files(expected_file, generated_file):
    """Compare BBH file formats."""
    print(f"\nComparing BBH file formats:")
    print(f"Expected: {expected_file}")
    print(f"Generated: {generated_file}")
    print("-" * 60)
    
    # Load files
    expected_df = pd.read_csv(expected_file)
    generated_df = pd.read_csv(generated_file)
    
    # Check columns
    expected_cols = list(expected_df.columns)
    generated_cols = list(generated_df.columns)
    
    print(f"\nExpected columns: {expected_cols}")
    print(f"Generated columns: {generated_cols}")
    
    if expected_cols == generated_cols:
        print("✓ Column names match!")
    else:
        print("✗ Column names do not match!")
        missing = set(expected_cols) - set(generated_cols)
        extra = set(generated_cols) - set(expected_cols)
        if missing:
            print(f"  Missing columns: {missing}")
        if extra:
            print(f"  Extra columns: {extra}")
    
    # Check data types
    print("\nData types comparison:")
    for col in expected_cols:
        if col in generated_df.columns:
            exp_dtype = expected_df[col].dtype
            gen_dtype = generated_df[col].dtype
            if exp_dtype == gen_dtype:
                print(f"  ✓ {col}: {exp_dtype}")
            else:
                print(f"  ✗ {col}: expected {exp_dtype}, got {gen_dtype}")
    
    # Check sample data
    print("\nSample data (first 3 rows):")
    print("\nExpected:")
    print(expected_df.head(3))
    print("\nGenerated:")
    print(generated_df.head(3))
    
    # Check BBH column values
    if 'BBH' in expected_df.columns and 'BBH' in generated_df.columns:
        expected_bbh_values = set(expected_df['BBH'].unique())
        generated_bbh_values = set(generated_df['BBH'].unique())
        print(f"\nBBH column values:")
        print(f"Expected: {expected_bbh_values}")
        print(f"Generated: {generated_bbh_values}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python test_bbh_format.py <expected_file> <generated_file>")
        sys.exit(1)
    
    compare_bbh_files(sys.argv[1], sys.argv[2])