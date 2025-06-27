"""
Multi-view ICA implementation for MultiModulon analysis.

This module provides GPU-accelerated multi-view ICA functionality
using PyTorch within Docker containers.
"""

import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd


def prepare_docker_environment(
    input_data_dir: str,
    output_dir: str,
    species_X_matrices: Dict[str, pd.DataFrame],
    a_values: Dict[str, int],
    c: int,
    mode: str = 'gpu'
) -> Tuple[str, str]:
    """
    Prepare the Docker environment by creating necessary files and directories.
    
    Args:
        input_data_dir: Directory containing input data
        output_dir: Directory for output files
        species_X_matrices: Dictionary mapping species names to X matrices
        a_values: Dictionary mapping species names to number of components
        c: Number of shared sources
        mode: 'gpu' or 'cpu' mode for computation
        
    Returns:
        Tuple of (temp_dir_path, script_path)
    """
    # Create temporary directory for Docker execution
    temp_dir = tempfile.mkdtemp(prefix="multiview_ica_")
    
    # Save X matrices as CSV files
    print("Saving X matrices for Docker execution...")
    data_dir = os.path.join(temp_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    
    species_list = list(species_X_matrices.keys())
    for species, X_matrix in species_X_matrices.items():
        csv_path = os.path.join(data_dir, f"{species}_X.csv")
        X_matrix.to_csv(csv_path, index=True)
    
    # Create the execution script
    script_content = f'''
import sys
import os
import pandas as pd
import numpy as np
import torch

# Add the multiview_ica module path
sys.path.append('/workspace')

# Import the functions from the original implementation
from multiview_ica_impl import run_multi_view_ICA_on_6_datasets

# Load the data
print("Loading X matrices...")
species_list = {species_list}
X_matrices = []
for species in species_list:
    csv_path = f"/workspace/data/{{species}}_X.csv"
    X = pd.read_csv(csv_path, index_col=0)
    X_matrices.append(X)
    print(f"Loaded {{species}}: {{X.shape}}")

# Get a_values in the same order
a_values = {[a_values[sp] for sp in species_list]}
c = {c}

print("\\nRunning multi-view ICA...")
print(f"a_values: {{a_values}}")
print(f"c (shared sources): {{c}}")

# Run the multi-view ICA
results = run_multi_view_ICA_on_6_datasets(
    *X_matrices,
    *a_values,
    c,
    mode='{mode}'
)

# Save results
print("\\nSaving results...")
for i, (species, result) in enumerate(zip(species_list, results)):
    output_path = f"/workspace/output/{{species}}_M.csv"
    result.to_csv(output_path, index=False)
    print(f"Saved {{species}} M matrix: {{result.shape}}")

print("\\nMulti-view ICA completed successfully!")
'''
    
    script_path = os.path.join(temp_dir, "run_multiview_ica.py")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Copy the original multiview_ica implementation
    impl_path = os.path.join(temp_dir, "multiview_ica_impl.py")
    with open("/home/gaoyuan/PhD/iMM_2/multimodulon_dev/multiview_ica.py", 'r') as src:
        with open(impl_path, 'w') as dst:
            dst.write(src.read())
    
    # Create output directory
    os.makedirs(os.path.join(temp_dir, "output"), exist_ok=True)
    
    return temp_dir, script_path


def run_multiview_ica_docker(
    species_X_matrices: Dict[str, pd.DataFrame],
    a_values: Dict[str, int],
    c: int,
    mode: str = 'gpu',
    docker_image: str = 'pytorch/pytorch:2.6.0-cuda12.4-cudnn9-runtime'
) -> Dict[str, pd.DataFrame]:
    """
    Run multi-view ICA using Docker with PyTorch.
    
    Args:
        species_X_matrices: Dictionary mapping species names to aligned X matrices
        a_values: Dictionary mapping species names to number of components
        c: Number of shared sources
        mode: 'gpu' or 'cpu' mode
        docker_image: Docker image to use
        
    Returns:
        Dictionary mapping species names to M matrices (ICA results)
    """
    print("\n" + "="*60)
    print("MULTI-VIEW ICA EXECUTION")
    print("="*60)
    
    start_time = time.time()
    
    # Validate inputs
    if len(species_X_matrices) != 6:
        raise ValueError(f"Expected exactly 6 species, got {len(species_X_matrices)}")
    
    if len(a_values) != 6:
        raise ValueError(f"Expected exactly 6 a_values, got {len(a_values)}")
    
    # Prepare Docker environment
    print("\n[1/4] Preparing Docker environment...")
    prep_start = time.time()
    temp_dir, script_path = prepare_docker_environment(
        input_data_dir="",
        output_dir="",
        species_X_matrices=species_X_matrices,
        a_values=a_values,
        c=c,
        mode=mode
    )
    print(f"✓ Environment prepared in {time.time() - prep_start:.2f}s")
    
    # Pull Docker image if needed
    print(f"\n[2/4] Checking Docker image: {docker_image}")
    pull_start = time.time()
    try:
        subprocess.run(
            ["docker", "image", "inspect", docker_image],
            capture_output=True,
            check=True
        )
        print(f"✓ Docker image already available")
    except subprocess.CalledProcessError:
        print(f"Pulling Docker image...")
        subprocess.run(
            ["docker", "pull", docker_image],
            check=True
        )
        print(f"✓ Docker image pulled in {time.time() - pull_start:.2f}s")
    
    # Run Docker container
    print(f"\n[3/4] Running multi-view ICA in Docker container...")
    run_start = time.time()
    
    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{temp_dir}:/workspace",
        "-w", "/workspace"
    ]
    
    # Add GPU support if in GPU mode
    if mode == 'gpu':
        docker_cmd.extend(["--gpus", "all"])
    
    docker_cmd.extend([
        docker_image,
        "python", "/workspace/run_multiview_ica.py"
    ])
    
    try:
        result = subprocess.run(
            docker_cmd,
            capture_output=True,
            text=True,
            check=True
        )
        print(result.stdout)
        if result.stderr:
            print("Warnings:", result.stderr)
        print(f"✓ Multi-view ICA completed in {time.time() - run_start:.2f}s")
    except subprocess.CalledProcessError as e:
        print(f"Error running Docker container:")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        raise
    
    # Load results
    print(f"\n[4/4] Loading results...")
    load_start = time.time()
    results = {}
    output_dir = os.path.join(temp_dir, "output")
    
    for species in species_X_matrices.keys():
        result_path = os.path.join(output_dir, f"{species}_M.csv")
        if os.path.exists(result_path):
            results[species] = pd.read_csv(result_path)
            print(f"✓ Loaded {species} M matrix: {results[species].shape}")
        else:
            raise FileNotFoundError(f"Result file not found: {result_path}")
    
    print(f"✓ Results loaded in {time.time() - load_start:.2f}s")
    
    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)
    
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"TOTAL EXECUTION TIME: {total_time:.2f}s")
    print(f"{'='*60}\n")
    
    return results


def run_multiview_ica_native(
    species_X_matrices: Dict[str, pd.DataFrame],
    a_values: Dict[str, int],
    c: int,
    mode: str = 'gpu'
) -> Dict[str, pd.DataFrame]:
    """
    Run multi-view ICA natively (without Docker).
    
    This requires PyTorch to be installed in the current environment.
    
    Args:
        species_X_matrices: Dictionary mapping species names to aligned X matrices
        a_values: Dictionary mapping species names to number of components
        c: Number of shared sources
        mode: 'gpu' or 'cpu' mode
        
    Returns:
        Dictionary mapping species names to M matrices (ICA results)
    """
    try:
        import torch
        sys.path.append('/home/gaoyuan/PhD/iMM_2/multimodulon_dev')
        from multiview_ica import run_multi_view_ICA_on_6_datasets
    except ImportError:
        raise ImportError(
            "PyTorch not found. Please use Docker mode or install PyTorch."
        )
    
    print("\n" + "="*60)
    print("MULTI-VIEW ICA EXECUTION (Native)")
    print("="*60)
    
    start_time = time.time()
    
    # Validate inputs
    if len(species_X_matrices) != 6:
        raise ValueError(f"Expected exactly 6 species, got {len(species_X_matrices)}")
    
    # Prepare data
    species_list = list(species_X_matrices.keys())
    X_matrices = [species_X_matrices[sp] for sp in species_list]
    a_list = [a_values[sp] for sp in species_list]
    
    print(f"\nRunning multi-view ICA...")
    print(f"Species: {species_list}")
    print(f"a_values: {a_list}")
    print(f"c (shared sources): {c}")
    print(f"Mode: {mode}")
    
    # Run multi-view ICA
    results = run_multi_view_ICA_on_6_datasets(
        *X_matrices,
        *a_list,
        c,
        mode=mode
    )
    
    # Create results dictionary
    results_dict = {}
    for species, result in zip(species_list, results):
        results_dict[species] = result
        print(f"✓ {species} M matrix: {result.shape}")
    
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"TOTAL EXECUTION TIME: {total_time:.2f}s")
    print(f"{'='*60}\n")
    
    return results_dict