"""
Multi-view ICA optimization functions for selecting the optimal number of core components.

This module implements the NRE (Noise Reduction Error) method for automatically
determining the optimal number of core components in multi-view ICA.
"""

import time
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Warning: matplotlib not available for plotting")




def calculate_nre_proper(
    S_matrices: List[pd.DataFrame],
    k_core: int
) -> float:
    """
    Calculate the Normalized Reconstruction Error (NRE) as defined in the paper.
    
    NRE(k) = (1/N_1) * Σᵢ₌₁^{N_1} Σ_{d=1}^D ||ẑ_d^{(k)}_i||² / k
    
    where:
    - ẑ_d^{(k)}_i = z_{d,0}^{(k)}_i - z̄_0^{(k)}_i (residual for sample i)
    - z_{d,0}^{(k)}_i are the first k components from view d for sample i  
    - z̄_0^{(k)}_i = (1/D) Σ_{ℓ=1}^D z_{ℓ,0}^{(k)}_i (cross-view mean for sample i)
    - N_1 is the number of test samples
    
    Args:
        S_matrices: List of source signal matrices with shape (samples, components)
        k_core: Number of core components to evaluate (k)
        
    Returns:
        NRE score (lower is better)
    """
    if k_core <= 0:
        return float('inf')
    
    D = len(S_matrices)  # Number of views
    N_1 = S_matrices[0].shape[0]  # Number of samples (rows)
    
    # Extract first k_core components from each view (z_{d,0}^(k))
    # S_matrices have shape (samples, components), so we need columns
    z_d_0_k = []
    for S in S_matrices:
        if k_core > S.shape[1]:  # Check number of components (columns)
            # Pad with zeros if k_core exceeds available components
            padded = np.zeros((S.shape[0], k_core))
            padded[:, :S.shape[1]] = S.values
            z_d_0_k.append(padded)
        else:
            z_d_0_k.append(S.values[:, :k_core])  # First k components (columns)
    
    # Calculate NRE for each sample, then average
    nre_per_sample = []
    
    for i in range(N_1):
        # Extract components for sample i from all views: z_{d,0}^{(k)}_i
        z_i_views = [z_d[i, :] for z_d in z_d_0_k]  # List of k_core arrays
        
        # Calculate z̄_0^{(k)}_i = (1/D) * Σ_{ℓ=1}^D z_{ℓ,0}^{(k)}_i
        z_bar_i = np.mean(z_i_views, axis=0)  # Shape: (k_core,)
        
        # Calculate Σ_{d=1}^D ||ẑ_d^{(k)}_i||² where ẑ_d^{(k)}_i = z_{d,0}^{(k)}_i - z̄_0^{(k)}_i
        sample_nre = 0.0
        for z_i_d in z_i_views:
            hat_z_i_d = z_i_d - z_bar_i  # Residual for this view
            sample_nre += np.sum(hat_z_i_d**2)  # ||ẑ_d^{(k)}_i||²
        
        # Divide by k for this sample: Σ_{d=1}^D ||ẑ_d^{(k)}_i||² / k
        nre_per_sample.append(sample_nre / k_core)
    
    # Average over all test samples: (1/N_1) * Σᵢ₌₁^{N_1} NRE(k)_i
    average_nre = np.mean(nre_per_sample)
    return average_nre






def run_nre_optimization(
    species_X_matrices: Dict[str, pd.DataFrame],
    k_candidates: List[int],
    max_a_per_view: int,
    train_frac: float = 0.75,
    num_runs: int = 3,
    mode: str = 'gpu',
    seed: int = 42
) -> Tuple[int, Dict[int, float], Dict[int, List], Dict[int, List], plt.Figure]:
    """
    NRE optimization using PyTorch.
    """
    
    species_list = list(species_X_matrices.keys())
    n_species = len(species_list)
    if n_species < 2:
        raise ValueError(f"Expected at least 2 species, got {n_species}")
    
    # Check and display device status
    import torch
    if mode == 'gpu' and torch.cuda.is_available():
        print(f"Using GPU: {torch.cuda.get_device_name(0)}")
    else:
        print("Using CPU")
    
    # Initialize storage for NRE results
    mean_nre_per_k = {}
    all_nre_per_k = {k: [] for k in k_candidates}
    
    for run in range(num_runs):
        
        # Split data consistently across views
        n_samples = species_X_matrices[species_list[0]].shape[0]  # samples are rows
        indices = np.arange(n_samples)
        train_idx, test_idx = train_test_split(
            indices, train_size=train_frac, random_state=seed + run
        )
        
        X_train_views = []
        X_test_views = []
        for species in species_list:
            X = species_X_matrices[species]
            X_train_views.append(X.iloc[train_idx, :])  # Select rows (samples)
            X_test_views.append(X.iloc[test_idx, :])    # Select rows (samples)
        
        for k in k_candidates:
            start_time = time.time()
            
            # Import the generalized function
            from .multiview_ica import run_multi_view_ICA_on_datasets
            
            # Run multi-view ICA with current k (training)
            train_sources = run_multi_view_ICA_on_datasets(
                X_train_views,
                [max_a_per_view] * n_species,
                k,
                mode=mode
            )
            
            # Train another model on test data to get comparable results
            test_sources = run_multi_view_ICA_on_datasets(
                X_test_views,
                [max_a_per_view] * n_species,
                k,
                mode=mode
            )
            
            # Calculate NRE using the test source signals (components x samples)
            nre = calculate_nre_proper(test_sources, k)
            all_nre_per_k[k].append(nre)
            
            elapsed_time = time.time() - start_time
            print(f"k={k}: {elapsed_time:.1f}s, NRE={nre:.6f}")
    
    # Calculate mean NRE for each k
    for k in k_candidates:
        if all_nre_per_k[k]:
            mean_nre_per_k[k] = np.mean(all_nre_per_k[k])
    
    # Find optimal k
    min_nre = min(mean_nre_per_k.values())
    tolerance = 1e-8
    optimal_k_candidates = [k for k, nre in mean_nre_per_k.items() 
                           if abs(nre - min_nre) < tolerance]
    best_k = max(optimal_k_candidates)
    
    print(f"\nOptimal k = {best_k} (NRE = {mean_nre_per_k[best_k]:.6f})")
    
    # Create NRE plot
    if 'plt' in globals():
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        # Create a dummy figure object if matplotlib is not available
        class DummyFig:
            def savefig(self, *args, **kwargs): pass
        fig = DummyFig()
        ax = None
    
    if ax is not None:
        k_values = sorted(mean_nre_per_k.keys())
        
        mean_values = [mean_nre_per_k[k] for k in k_values]
        std_values = [np.std(all_nre_per_k[k]) if len(all_nre_per_k[k]) > 1 else 0 for k in k_values]
        
        # Plot NRE with error bars
        ax.errorbar(k_values, mean_values, yerr=std_values, 
                   marker='o', markersize=8, capsize=5, capthick=2,
                   color='blue', label='NRE Score', linewidth=2)
        
        # Highlight optimal k
        ax.axvline(x=best_k, color='red', linestyle='--', alpha=0.7, linewidth=2,
                  label=f'Optimal k = {best_k}')
        ax.scatter([best_k], [mean_nre_per_k[best_k]], color='red', s=120, zorder=5,
                  edgecolors='darkred', linewidth=2)
        
        ax.set_xlabel('Number of Core Components (k)', fontsize=12)
        ax.set_ylabel('NRE Score (Lower is Better)', fontsize=12)
        ax.set_title('Normalized Reconstruction Error (NRE) vs Number of Core Components', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add optimal value text
        ax.text(0.02, 0.98, f'Optimal k = {best_k}\nNRE = {mean_nre_per_k[best_k]:.6f}', 
               transform=ax.transAxes, fontsize=11, fontweight='bold',
               bbox=dict(boxstyle="round,pad=0.4", facecolor='lightblue', alpha=0.8),
               verticalalignment='top')
        
        # Adjust layout
        plt.tight_layout()
    
    return best_k, mean_nre_per_k, all_nre_per_k, fig