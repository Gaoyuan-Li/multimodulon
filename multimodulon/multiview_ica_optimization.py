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




def calculate_rgr(
    S_matrices: List[pd.DataFrame],
    k_core: int,
    k_prev: Optional[int] = None,
    prev_reconstruction_error: Optional[float] = None
) -> float:
    """
    Calculate the Reconstruction Gain Rate (RGR).
    
    RGR measures the relative improvement in reconstruction quality per added component:
    RGR(k) = (RE(k-Δk) - RE(k)) / (Δk * RE(k))
    
    where:
    - RE(k) is the reconstruction error with k components
    - Δk is the step size (default 5)
    - Higher RGR indicates better gain from adding components
    
    Args:
        S_matrices: List of source signal matrices with shape (samples, components)
        k_core: Number of core components to evaluate
        k_prev: Previous k value (for calculating gain)
        prev_reconstruction_error: RE at k_prev (to avoid recalculation)
        
    Returns:
        RGR score (higher is better, optimal at elbow point)
    """
    # Calculate reconstruction error for current k
    re_k = calculate_reconstruction_error(S_matrices, k_core)
    
    if k_prev is None or prev_reconstruction_error is None:
        # For first k, return normalized RE inverse as proxy
        return 1.0 / (re_k + 1e-8)
    
    # Calculate gain rate
    delta_k = k_core - k_prev
    if delta_k <= 0 or re_k >= prev_reconstruction_error:
        return 0.0  # No improvement
    
    rgr = (prev_reconstruction_error - re_k) / (delta_k * re_k)
    return rgr


def calculate_reconstruction_error(
    S_matrices: List[pd.DataFrame],
    k_core: int
) -> float:
    """
    Calculate the reconstruction error using first k components.
    
    RE(k) = Σᵢ Σ_d ||X_d_i - X̂_d_i||² / (N * D)
    
    where X̂_d_i is reconstructed from k core components
    """
    if k_core <= 0:
        return float('inf')
    
    D = len(S_matrices)
    N = S_matrices[0].shape[0]
    
    total_error = 0.0
    
    # Extract first k components from each view
    core_components = []
    for S in S_matrices:
        if k_core > S.shape[1]:
            padded = np.zeros((S.shape[0], k_core))
            padded[:, :S.shape[1]] = S.values
            core_components.append(padded)
        else:
            core_components.append(S.values[:, :k_core])
    
    # Calculate reconstruction error
    for d, S in enumerate(S_matrices):
        # Use full matrix for error calculation
        full_signal = S.values
        core_signal = core_components[d]
        
        # Simple reconstruction error: difference from using only k components
        # In practice, this measures how much signal is in non-core components
        if k_core < S.shape[1]:
            non_core_signal = S.values[:, k_core:]
            error = np.sum(non_core_signal ** 2)
        else:
            error = 0.0
        
        total_error += error
    
    return total_error / (N * D)


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
    seed: int = 42,
    metric: str = 'nre'
) -> Tuple[int, Dict[int, float], Dict[int, List], Dict[int, List], plt.Figure]:
    """
    Optimization using NRE or RGR metric.
    
    Args:
        metric: 'nre' for Normalized Reconstruction Error, 'rgr' for Reconstruction Gain Rate
    """
    
    species_list = list(species_X_matrices.keys())
    n_species = len(species_list)
    if n_species < 2:
        raise ValueError(f"Expected at least 2 species, got {n_species}")
    
    if metric not in ['nre', 'rgr']:
        raise ValueError(f"Unknown metric: {metric}. Use 'nre' or 'rgr'")
    
    # Check and display device status
    import torch
    if mode == 'gpu' and torch.cuda.is_available():
        print(f"Using GPU: {torch.cuda.get_device_name(0)}")
    else:
        print("Using CPU")
    
    # Initialize storage for metric results
    mean_metric_per_k = {}
    all_metric_per_k = {k: [] for k in k_candidates}
    
    # For RGR, we need to track previous reconstruction errors
    if metric == 'rgr':
        prev_re_per_run = {run: {} for run in range(num_runs)}
    
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
            
            # Calculate metric using the test source signals (components x samples)
            if metric == 'nre':
                score = calculate_nre_proper(test_sources, k)
            else:  # rgr
                # Find previous k value for RGR calculation
                k_idx = k_candidates.index(k)
                if k_idx > 0:
                    k_prev = k_candidates[k_idx - 1]
                    prev_re = prev_re_per_run[run].get(k_prev)
                    score = calculate_rgr(test_sources, k, k_prev, prev_re)
                else:
                    score = calculate_rgr(test_sources, k)
                
                # Store reconstruction error for next iteration
                prev_re_per_run[run][k] = calculate_reconstruction_error(test_sources, k)
            
            all_metric_per_k[k].append(score)
            
            elapsed_time = time.time() - start_time
            metric_name = metric.upper()
            print(f"k={k}: {elapsed_time:.1f}s, {metric_name}={score:.6f}")
    
    # Calculate mean metric for each k
    for k in k_candidates:
        if all_metric_per_k[k]:
            mean_metric_per_k[k] = np.mean(all_metric_per_k[k])
    
    # Find optimal k
    if metric == 'nre':
        # For NRE, minimize
        optimal_value = min(mean_metric_per_k.values())
        tolerance = 1e-8
        optimal_k_candidates = [k for k, score in mean_metric_per_k.items() 
                               if abs(score - optimal_value) < tolerance]
        best_k = max(optimal_k_candidates)  # Choose largest k if tied
    else:  # rgr
        # For RGR, find elbow point (maximum gain rate)
        optimal_value = max(mean_metric_per_k.values())
        tolerance = optimal_value * 0.1  # 10% tolerance for RGR
        optimal_k_candidates = [k for k, score in mean_metric_per_k.items() 
                               if score >= optimal_value - tolerance]
        best_k = min(optimal_k_candidates)  # Choose smallest k if tied
    
    metric_name = metric.upper()
    print(f"\nOptimal k = {best_k} ({metric_name} = {mean_metric_per_k[best_k]:.6f})")
    
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
        k_values = sorted(mean_metric_per_k.keys())
        
        mean_values = [mean_metric_per_k[k] for k in k_values]
        std_values = [np.std(all_metric_per_k[k]) if len(all_metric_per_k[k]) > 1 else 0 for k in k_values]
        
        # Plot metric with error bars
        metric_label = 'NRE Score' if metric == 'nre' else 'RGR Score'
        ax.errorbar(k_values, mean_values, yerr=std_values, 
                   marker='o', markersize=8, capsize=5, capthick=2,
                   color='blue', label=metric_label, linewidth=2)
        
        # Highlight optimal k
        ax.axvline(x=best_k, color='red', linestyle='--', alpha=0.7, linewidth=2,
                  label=f'Optimal k = {best_k}')
        ax.scatter([best_k], [mean_metric_per_k[best_k]], color='red', s=120, zorder=5,
                  edgecolors='darkred', linewidth=2)
        
        ax.set_xlabel('Number of Core Components (k)', fontsize=12)
        if metric == 'nre':
            ax.set_ylabel('NRE Score (Lower is Better)', fontsize=12)
            ax.set_title('Normalized Reconstruction Error (NRE) vs Number of Core Components', fontsize=14, fontweight='bold')
        else:
            ax.set_ylabel('RGR Score (Higher is Better)', fontsize=12)
            ax.set_title('Reconstruction Gain Rate (RGR) vs Number of Core Components', fontsize=14, fontweight='bold')
        
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add optimal value text
        metric_name = metric.upper()
        ax.text(0.02, 0.98, f'Optimal k = {best_k}\n{metric_name} = {mean_metric_per_k[best_k]:.6f}', 
               transform=ax.transAxes, fontsize=11, fontweight='bold',
               bbox=dict(boxstyle="round,pad=0.4", facecolor='lightblue', alpha=0.8),
               verticalalignment='top')
        
        # Adjust layout
        plt.tight_layout()
    
    return best_k, mean_metric_per_k, all_metric_per_k, fig