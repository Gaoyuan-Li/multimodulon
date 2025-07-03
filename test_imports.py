#!/usr/bin/env python3
"""Test script to verify all imports work correctly after reorganization."""

print("Testing imports...")

# Test basic import
try:
    from multimodulon import MultiModulon
    print("✓ Basic import works: from multimodulon import MultiModulon")
except ImportError as e:
    print(f"✗ Basic import failed: {e}")

# Test individual module imports
modules = [
    "multimodulon.core",
    "multimodulon.core_io", 
    "multimodulon.plotting",
    "multimodulon.gene_alignment",
    "multimodulon.multiview_ica",
    "multimodulon.optimization",
    "multimodulon.utils",
    "multimodulon.species_data",
    "multimodulon.gff_utils"
]

for module in modules:
    try:
        __import__(module)
        print(f"✓ Import works: {module}")
    except ImportError as e:
        print(f"✗ Import failed: {module} - {e}")

# Test specific function imports
try:
    from multimodulon.utils import BBHAnalyzer, gff2pandas, extract_protein_sequences
    print("✓ Utils functions import correctly")
except ImportError as e:
    print(f"✗ Utils functions import failed: {e}")

try:
    from multimodulon.core_io import save_to_json_multimodulon, load_json_multimodulon
    print("✓ IO functions import correctly")
except ImportError as e:
    print(f"✗ IO functions import failed: {e}")

try:
    from multimodulon.plotting import view_iModulon_weights, compare_core_iModulon
    print("✓ Plotting functions import correctly")
except ImportError as e:
    print(f"✗ Plotting functions import failed: {e}")

print("\nImport test complete!")