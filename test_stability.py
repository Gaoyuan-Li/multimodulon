#!/usr/bin/env python
"""Test script for core_iModulon_stability function."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from multimodulon import MultiModulon

# Test with sample data
def test_stability():
    print("Testing core_iModulon_stability function...")
    print("=" * 60)
    
    # Initialize MultiModulon (assuming you have data loaded)
    # Update this path to your actual data
    input_path = "../path_to_your_data/Input_Data"  
    
    try:
        # Load multimodulon
        print("Loading MultiModulon...")
        mm = MultiModulon(input_path)
        
        # Check if M matrices exist
        has_m_matrix = False
        for species in mm.species:
            if mm[species]._M is not None:
                has_m_matrix = True
                break
        
        if not has_m_matrix:
            print("No M matrices found. Please run ICA first.")
            return
        
        # Find a core component to test
        test_component = None
        for species in mm.species:
            if mm[species]._M is not None:
                for col in mm[species].M.columns:
                    if col.startswith('Core_'):
                        test_component = col
                        break
                if test_component:
                    break
        
        if not test_component:
            print("No Core components found in M matrices.")
            return
        
        print(f"\nTesting with component: {test_component}")
        print("-" * 40)
        
        # Test different threshold methods
        methods = ["mad", "clustering", "otsu", "elbow"]
        
        for method in methods:
            print(f"\nTesting {method} method:")
            try:
                stable, threshold, scores = mm.core_iModulon_stability(
                    test_component,
                    threshold_method=method,
                    show_stats=True
                )
                
                print(f"  Threshold: {threshold:.3f}")
                print(f"  Stable species: {stable}")
                print(f"  Scores: {scores}")
                
            except Exception as e:
                print(f"  Error: {e}")
        
        # Test manual threshold
        print("\nTesting manual threshold (0.7):")
        try:
            stable, threshold, scores = mm.core_iModulon_stability(
                test_component,
                threshold_method="manual",
                manual_threshold=0.7,
                show_stats=True
            )
            
            print(f"  Threshold: {threshold:.3f}")
            print(f"  Stable species: {stable}")
            
        except Exception as e:
            print(f"  Error: {e}")
        
        print("\n" + "=" * 60)
        print("Testing completed successfully!")
        
    except Exception as e:
        print(f"Error loading MultiModulon: {e}")
        print("\nPlease update the input_path variable to point to your data.")


if __name__ == "__main__":
    # Create a simple test
    print("Creating a simple test for the stability function...")
    print("\nTo run a full test, please:")
    print("1. Update the input_path in this script to point to your data")
    print("2. Ensure you have M matrices computed (run ICA first)")
    print("3. Run this script again")
    
    # You can uncomment this to run the test with your data
    # test_stability()