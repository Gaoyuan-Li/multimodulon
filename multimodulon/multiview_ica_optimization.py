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
from tqdm import tqdm

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Warning: matplotlib not available for plotting")

# Cohen's d effect size calculation for component selection




def calculate_cohens_d_effect_size(weight_vector: np.ndarray, seed: int = 42, num_top_gene: int = 20) -> float:
    """
    Calculate Cohen's d effect size between top N genes and the rest.
    
    Cohen's d = (mean_topN - mean_rest) / pooled_std
    
    Args:
        weight_vector: 1D array of weights for a single component
        seed: Random seed (kept for compatibility but not used)
        num_top_gene: Number of top genes to use for effect size calculation
        
    Returns:
        Cohen's d effect size measuring separation between top N genes and rest
    """
    if len(weight_vector) < num_top_gene + 1:  # Need at least num_top_gene + 1 genes
        return 0.0
    
    # Get absolute values and indices of top N genes
    abs_weights = np.abs(weight_vector)
    top_n_indices = np.argpartition(abs_weights, -num_top_gene)[-num_top_gene:]
    rest_indices = np.setdiff1d(np.arange(len(weight_vector)), top_n_indices)
    
    # Get the actual weight values for each group
    top_n_weights = weight_vector[top_n_indices]
    rest_weights = weight_vector[rest_indices]
    
    # Calculate means
    mean_top_n = np.mean(top_n_weights)
    mean_rest = np.mean(rest_weights)
    
    # Calculate standard deviations
    std_top_n = np.std(top_n_weights, ddof=1)
    std_rest = np.std(rest_weights, ddof=1)
    
    # Sample sizes
    n1 = len(top_n_weights)
    n2 = len(rest_weights)
    
    # Calculate pooled standard deviation
    pooled_std = np.sqrt(((n1 - 1) * std_top_n**2 + (n2 - 1) * std_rest**2) / (n1 + n2 - 2))
    
    # Avoid division by zero
    if pooled_std == 0:
        return 0.0
    
    # Calculate Cohen's d
    cohens_d = abs(mean_top_n - mean_rest) / pooled_std
    
    return cohens_d


def calculate_average_effect_sizes(M_matrices: Dict[str, pd.DataFrame], seed: int = 42, num_top_gene: int = 20) -> List[float]:
    """
    Calculate mean effect sizes across species for each component.
    
    Args:
        M_matrices: Dict mapping species names to M matrices (samples, components)
        num_top_gene: Number of top genes to use for effect size calculation
        
    Returns:
        List of mean effect sizes, one per component
    """
    # Sort species list for consistent ordering across runs
    species_list = sorted(list(M_matrices.keys()))
    n_components = M_matrices[species_list[0]].shape[1]
    
    avg_effect_sizes = []
    
    for comp_idx in range(n_components):
        component_effect_sizes = []
        
        # Calculate effect size for this component across all species
        for species in species_list:
            M_matrix = M_matrices[species]
            weight_vector = M_matrix.iloc[:, comp_idx].values
            effect_size = calculate_cohens_d_effect_size(weight_vector, seed, num_top_gene)
            component_effect_sizes.append(effect_size)
        
        # Mean across species
        mean_effect_size = np.mean(component_effect_sizes)
        avg_effect_sizes.append(mean_effect_size)
    
    return avg_effect_sizes


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
    num_runs: int = 1,
    mode: str = 'gpu',
    seed: int = 42,
    metric: str = 'nre',
    threshold: Optional[float] = None,
    effective_size_threshold: float = 5,
    num_top_gene: int = 20
) -> Tuple[int, Dict[int, float], Dict[int, List], Dict[int, List], plt.Figure]:
    """
    Optimization using NRE or Cohen's d metric.
    
    Args:
        metric: 'nre' for Normalized Reconstruction Error, 'effect_size' for Cohen's d effect size
    """
    
    # Sort species list for consistent ordering across runs
    species_list = sorted(list(species_X_matrices.keys()))
    n_species = len(species_list)
    if n_species < 2:
        raise ValueError(f"Expected at least 2 species, got {n_species}")
    
    if metric not in ['nre', 'effect_size']:
        raise ValueError(f"Unknown metric: {metric}. Use 'nre' or 'effect_size' (Cohen's d)")
    
    # Set all random seeds for reproducibility (crucial after kernel restart)
    import torch
    import random
    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        # Additional CUDA determinism
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    
    # Check and display device status
    if mode == 'gpu' and torch.cuda.is_available():
        print(f"Using GPU: {torch.cuda.get_device_name(0)}")
    else:
        print("Using CPU")
    
    # Initialize storage for metric results
    mean_metric_per_k = {}
    all_metric_per_k = {k: [] for k in k_candidates}
    
    # For Cohen's d, also store individual component effect sizes for threshold analysis
    component_effect_sizes_per_k = {k: [] for k in k_candidates} if metric == 'effect_size' else None
    # For the new logic, store number of components above threshold
    num_above_threshold_per_k = {} if metric == 'effect_size' else None
    
    for run in range(num_runs):
        
        if metric == 'nre':
            # Split data consistently across views for NRE
            n_samples = species_X_matrices[species_list[0]].shape[0]  # samples are rows
            indices = np.arange(n_samples)
            train_idx, test_idx = train_test_split(
                indices, train_size=train_frac, random_state=seed
            )
            
            X_train_views = []
            X_test_views = []
            for species in species_list:
                X = species_X_matrices[species]
                X_train_views.append(X.iloc[train_idx, :])  # Select rows (samples)
                X_test_views.append(X.iloc[test_idx, :])    # Select rows (samples)
        else:
            # For effect_size, use full dataset (no train/test split)
            X_train_views = None  # Not used for effect_size
            X_test_views = None   # Not used for effect_size
        
        # Use tqdm for progress bar
        k_iterator = tqdm(k_candidates, desc=f"Run {run+1}/{num_runs}") if metric == 'effect_size' else k_candidates
        
        for k in k_iterator:
            start_time = time.time()
            
            # Import the generalized function
            from .multiview_ica import run_multi_view_ICA_on_datasets
            
            if metric == 'nre':
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
                
                # Calculate NRE using the test source signals
                score = calculate_nre_proper(test_sources, k)
                
            else:  # effect_size
                # For Cohen's d, we need the M matrices (unmixing matrices)
                from .multiview_ica import run_multiview_ica
                
                # For effect_size metric, use square ICA (a=c=k)
                a_values = {species: k for species in species_list}
                
                # Use full dataset for effect_size calculation
                M_matrices = run_multiview_ica(
                    species_X_matrices,  # Use the full dataset directly
                    a_values,
                    k,  # c = k (number of core components, square ICA)
                    mode=mode,
                    seed=seed  # Use exact same seed for reproducibility
                )
                
                # Calculate mean effect sizes ONLY for the first k (core) components
                core_avg_effect_sizes = []
                for comp_idx in range(k):  # Only first k components are core
                    component_effect_sizes = []
                    for species in species_list:
                        M_matrix = M_matrices[species]
                        weight_vector = M_matrix.iloc[:, comp_idx].values
                        effect_size = calculate_cohens_d_effect_size(weight_vector, seed, num_top_gene)
                        component_effect_sizes.append(effect_size)
                    # Mean across species
                    mean_effect_size = np.mean(component_effect_sizes)
                    core_avg_effect_sizes.append(mean_effect_size)
                
                # Store individual component effect sizes for threshold analysis
                if component_effect_sizes_per_k is not None:
                    component_effect_sizes_per_k[k].append(core_avg_effect_sizes)
                
                # Count CORE components above effective_size_threshold
                num_above_threshold = sum(1 for effect_size in core_avg_effect_sizes if effect_size > effective_size_threshold)
                
                # Store the count for this k
                if num_above_threshold_per_k is not None:
                    if k not in num_above_threshold_per_k:
                        num_above_threshold_per_k[k] = []
                    num_above_threshold_per_k[k].append(num_above_threshold)
                
                # Use number of components above threshold as the metric
                score = num_above_threshold
            
            all_metric_per_k[k].append(score)
            
            elapsed_time = time.time() - start_time
            if metric == 'nre':
                print(f"k={k}: time={elapsed_time:.1f}s, NRE={score:.6f}")
    
    # Calculate mean metric for each k
    for k in k_candidates:
        if all_metric_per_k[k]:
            # For both metrics, take mean of runs
            mean_metric_per_k[k] = np.mean(all_metric_per_k[k])
    
    # Find optimal k
    if metric == 'nre':
        # For NRE, minimize
        optimal_value = min(mean_metric_per_k.values())
        tolerance = 1e-8
        optimal_k_candidates = [k for k, score in mean_metric_per_k.items() 
                               if abs(score - optimal_value) < tolerance]
        best_k = max(optimal_k_candidates)
    else:  # effect_size
        # For effect_size, find where the growth plateaus
        k_sorted = sorted(mean_metric_per_k.keys())
        
        # Find the k where growth stops or plateaus
        best_k = k_sorted[0]  # Default to smallest k
        
        for i in range(len(k_sorted) - 1):
            current_k = k_sorted[i]
            next_k = k_sorted[i + 1]
            
            current_count = mean_metric_per_k[current_k]
            next_count = mean_metric_per_k[next_k]
            
            # If the count doesn't increase for the next k, we've found our plateau
            if next_count <= current_count:
                best_k = current_k
                break
            # If we reach the end and it's still growing, take the last k
            elif i == len(k_sorted) - 2:
                best_k = next_k
    
    if metric == 'nre':
        print(f"\nOptimal k = {best_k} (NRE = {mean_metric_per_k[best_k]:.6f})")
    else:
        print(f"\nOptimal k = {best_k} (Number of components above threshold = {mean_metric_per_k[best_k]:.1f})")
    
    # Threshold analysis for effect_size metric
    if metric == 'effect_size' and threshold is not None and component_effect_sizes_per_k is not None:
        print(f"\nThreshold Analysis (threshold = {threshold:.3f}):")
        print("-" * 50)
        
        for k in sorted(k_candidates):
            if component_effect_sizes_per_k[k]:
                # Average effect sizes across runs for each component
                all_runs_effects = component_effect_sizes_per_k[k]
                n_components = len(all_runs_effects[0])
                
                # Average each component across runs
                avg_component_effects = []
                for comp_idx in range(n_components):
                    comp_effects_across_runs = [run_effects[comp_idx] for run_effects in all_runs_effects]
                    avg_component_effects.append(np.mean(comp_effects_across_runs))
                
                # Count components above threshold
                above_threshold = [eff for eff in avg_component_effects if eff > threshold]
                n_above = len(above_threshold)
                percentage = (n_above / n_components) * 100
                
                print(f"k={k:2d}: {n_above:2d}/{n_components:2d} components above threshold ({percentage:5.1f}%)")
        
        print("-" * 50)
    
    # Create plot
    if 'plt' in globals():
        # Create line plot for both metrics
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        # Create a dummy figure object if matplotlib is not available
        class DummyFig:
            def savefig(self, *args, **kwargs): pass
        fig = DummyFig()
        ax = None
    
    if ax is not None:
        k_values = sorted(mean_metric_per_k.keys())
        
        if metric == 'nre':
            # Standard line plot for NRE
            mean_values = [mean_metric_per_k[k] for k in k_values]
            std_values = [np.std(all_metric_per_k[k]) if len(all_metric_per_k[k]) > 1 else 0 for k in k_values]
            
            ax.errorbar(k_values, mean_values, yerr=std_values, 
                       marker='o', markersize=8, capsize=5, capthick=2,
                       color='blue', label='NRE Score', linewidth=2)
            
            # Highlight optimal k
            ax.axvline(x=best_k, color='red', linestyle='--', alpha=0.7, linewidth=2,
                      label=f'Optimal k = {best_k}')
            ax.scatter([best_k], [mean_metric_per_k[best_k]], color='red', s=120, zorder=5,
                      edgecolors='darkred', linewidth=2)
            
            ax.set_ylabel('NRE Score (Lower is Better)', fontsize=12)
            ax.set_title('Normalized Reconstruction Error (NRE) vs Number of Core Components', fontsize=14, fontweight='bold')
            
        else:
            # Line plot for effect_size showing number of components above threshold
            mean_values = [mean_metric_per_k[k] for k in k_values]
            
            ax.plot(k_values, mean_values, 'go-', linewidth=2, markersize=8, 
                   label='Number of components above threshold')
            
            # Highlight optimal k
            ax.axvline(x=best_k, color='red', linestyle='--', alpha=0.7, linewidth=2,
                      label=f'Optimal k = {best_k}')
            ax.scatter([best_k], [mean_metric_per_k[best_k]], color='red', s=120, zorder=5,
                      edgecolors='darkred', linewidth=2)
            
            ax.set_ylabel('Number of Components Above Threshold', fontsize=12)
            ax.set_title(f'Number of Components with Cohen\'s d > {effective_size_threshold} vs Number of Core Components', 
                        fontsize=14, fontweight='bold')
        
        ax.set_xlabel('Number of Core Components (k)', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add optimal value text
        if metric == 'nre':
            label_text = f'Optimal k = {best_k}\nNRE = {mean_metric_per_k[best_k]:.6f}'
        else:
            label_text = f'Optimal k = {best_k}\n{int(mean_metric_per_k[best_k])} components above threshold'
        
        ax.text(0.02, 0.98, label_text, 
               transform=ax.transAxes, fontsize=11, fontweight='bold',
               bbox=dict(boxstyle="round,pad=0.4", facecolor='lightblue', alpha=0.8),
               verticalalignment='top')
        
        # Adjust layout
        plt.tight_layout()
    
    return best_k, mean_metric_per_k, all_metric_per_k, fig