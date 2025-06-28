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

try:
    from sklearn.mixture import GaussianMixture
except ImportError:
    print("Warning: sklearn not available for GMM analysis")




def calculate_gmm_effect_size(weight_vector: np.ndarray) -> float:
    """
    Calculate effect size from 2-component GMM fit to weight vector.
    
    Effect Size = |μ₁ - μ₂| / max(σ₁, σ₂)
    
    Args:
        weight_vector: 1D array of weights for a single component
        
    Returns:
        Effect size measuring distinction between central and tail groups
    """
    if len(weight_vector) < 4:  # Need minimum samples for GMM
        return 0.0
    
    # Reshape for sklearn
    X = weight_vector.reshape(-1, 1)
    
    try:
        # Fit 2-component GMM with fixed random state for reproducibility
        gmm = GaussianMixture(n_components=2, random_state=42, max_iter=100, init_params='kmeans')
        gmm.fit(X)
        
        # Extract means and standard deviations
        means = gmm.means_.flatten()
        covs = gmm.covariances_.flatten()
        stds = np.sqrt(covs)
        
        # Calculate effect size
        mu1, mu2 = means[0], means[1]
        sigma1, sigma2 = stds[0], stds[1]
        
        effect_size = abs(mu1 - mu2) / max(sigma1, sigma2)
        
        return effect_size
        
    except Exception:
        # GMM failed to converge or other issues
        return 0.0


def calculate_average_effect_sizes(M_matrices: Dict[str, pd.DataFrame]) -> List[float]:
    """
    Calculate mean effect sizes across species for each component.
    
    Args:
        M_matrices: Dict mapping species names to M matrices (samples, components)
        
    Returns:
        List of mean effect sizes, one per component
    """
    species_list = list(M_matrices.keys())
    n_components = M_matrices[species_list[0]].shape[1]
    
    avg_effect_sizes = []
    
    for comp_idx in range(n_components):
        component_effect_sizes = []
        
        # Calculate effect size for this component across all species
        for species in species_list:
            M_matrix = M_matrices[species]
            weight_vector = M_matrix.iloc[:, comp_idx].values
            effect_size = calculate_gmm_effect_size(weight_vector)
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
    threshold: Optional[float] = None
) -> Tuple[int, Dict[int, float], Dict[int, List], Dict[int, List], plt.Figure]:
    """
    Optimization using NRE or GMM metric.
    
    Args:
        metric: 'nre' for Normalized Reconstruction Error, 'gmm' for GMM effect size
    """
    
    # Sort species list for consistent ordering across runs
    species_list = sorted(list(species_X_matrices.keys()))
    n_species = len(species_list)
    if n_species < 2:
        raise ValueError(f"Expected at least 2 species, got {n_species}")
    
    if metric not in ['nre', 'gmm']:
        raise ValueError(f"Unknown metric: {metric}. Use 'nre' or 'gmm'")
    
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
    
    # For GMM, also store individual component effect sizes for threshold analysis
    component_effect_sizes_per_k = {k: [] for k in k_candidates} if metric == 'gmm' else None
    
    for run in range(num_runs):
        
        # Split data consistently across views
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
        
        for k in k_candidates:
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
                
            else:  # gmm
                # For GMM, we need the M matrices (unmixing matrices)
                # Run multi-view ICA with a=c=k (square ICA)
                from .multiview_ica import run_multiview_ica
                
                # Prepare a_values dict for GMM (a=c=k for all species)
                a_values = {species: k for species in species_list}
                
                # Run on test data to get M matrices
                M_matrices = run_multiview_ica(
                    {species: X_test_views[i] for i, species in enumerate(species_list)},
                    a_values,
                    k,  # c = k
                    mode=mode,
                    seed=seed  # Use exact same seed for reproducibility
                )
                
                # Calculate mean effect sizes across components
                avg_effect_sizes = calculate_average_effect_sizes(M_matrices)
                
                # Store individual component effect sizes for threshold analysis
                if component_effect_sizes_per_k is not None:
                    component_effect_sizes_per_k[k].append(avg_effect_sizes)
                
                # Use median of all component effect sizes as the metric
                score = np.median(avg_effect_sizes)
            
            all_metric_per_k[k].append(score)
            
            elapsed_time = time.time() - start_time
            if metric == 'nre':
                print(f"k={k}: time={elapsed_time:.1f}s, NRE={score:.6f}")
            else:  # gmm
                print(f"k={k}: time={elapsed_time:.1f}s, median GMM effect size={score:.6f}")
    
    # Calculate mean metric for each k (for NRE) or median metric (for GMM)
    for k in k_candidates:
        if all_metric_per_k[k]:
            if metric == 'gmm':
                # For GMM, we already stored median values, just take mean of runs
                mean_metric_per_k[k] = np.mean(all_metric_per_k[k])
            else:
                # For NRE, take mean as before
                mean_metric_per_k[k] = np.mean(all_metric_per_k[k])
    
    # Find optimal k
    if metric == 'nre':
        # For NRE, minimize
        optimal_value = min(mean_metric_per_k.values())
        tolerance = 1e-8
        optimal_k_candidates = [k for k, score in mean_metric_per_k.items() 
                               if abs(score - optimal_value) < tolerance]
        best_k = max(optimal_k_candidates)
    else:  # gmm
        # For GMM effect size, maximize
        optimal_value = max(mean_metric_per_k.values())
        tolerance = optimal_value * 0.05  # 5% tolerance
        optimal_k_candidates = [k for k, score in mean_metric_per_k.items() 
                               if score >= optimal_value - tolerance]
        best_k = min(optimal_k_candidates)  # Choose smallest k if tied
    
    if metric == 'nre':
        print(f"\nOptimal k = {best_k} (NRE = {mean_metric_per_k[best_k]:.6f})")
    else:
        print(f"\nOptimal k = {best_k} (Median GMM effect size = {mean_metric_per_k[best_k]:.6f})")
    
    # Threshold analysis for GMM metric
    if metric == 'gmm' and threshold is not None and component_effect_sizes_per_k is not None:
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
        if metric == 'gmm' and component_effect_sizes_per_k is not None:
            # Create box plot for GMM showing distribution of component effect sizes
            fig, ax = plt.subplots(figsize=(14, 8))
        else:
            # Regular line plot for NRE
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
            # Box plot for GMM showing distribution of component effect sizes
            if component_effect_sizes_per_k is not None:
                # Prepare data for box plot
                box_data = []
                box_positions = []
                
                for k in k_values:
                    if component_effect_sizes_per_k[k]:
                        # Get all component effect sizes for this k (averaged across runs)
                        all_runs_effects = component_effect_sizes_per_k[k]
                        n_components = len(all_runs_effects[0])
                        
                        # Average each component across runs
                        avg_component_effects = []
                        for comp_idx in range(n_components):
                            comp_effects_across_runs = [run_effects[comp_idx] for run_effects in all_runs_effects]
                            avg_component_effects.append(np.mean(comp_effects_across_runs))
                        
                        box_data.append(avg_component_effects)
                        box_positions.append(k)
                
                # Create box plot
                bp = ax.boxplot(box_data, positions=box_positions, widths=2.5, patch_artist=True,
                               boxprops=dict(facecolor='lightgreen', alpha=0.7),
                               medianprops=dict(color='darkgreen', linewidth=2),
                               whiskerprops=dict(color='darkgreen'),
                               capprops=dict(color='darkgreen'),
                               flierprops=dict(marker='o', markerfacecolor='red', markersize=4, alpha=0.5))
                
                # Add mean line
                mean_values = [mean_metric_per_k[k] for k in k_values]
                ax.plot(k_values, mean_values, 'o-', linewidth=2, markersize=6, 
                       label='Mean Effect Size', color='red', alpha=0.8)
                
                # Add threshold line if specified
                if threshold is not None:
                    ax.axhline(y=threshold, color='orange', linestyle=':', alpha=0.8, linewidth=2,
                              label=f'Threshold = {threshold:.3f}')
                
                ax.set_ylabel('GMM Effect Size (Higher is Better)', fontsize=12)
                ax.set_title('Distribution of GMM Effect Sizes per Component vs Number of Core Components', fontsize=14, fontweight='bold')
                
                # Highlight optimal k
                ax.axvline(x=best_k, color='red', linestyle='--', alpha=0.7, linewidth=2,
                          label=f'Optimal k = {best_k}')
            
            else:
                # Fallback to line plot if no component data available
                mean_values = [mean_metric_per_k[k] for k in k_values]
                ax.plot(k_values, mean_values, 'go-', linewidth=2, markersize=6, label='Mean GMM Effect Size')
                ax.set_ylabel('Mean GMM Effect Size (Higher is Better)', fontsize=12)
                ax.set_title('Mean GMM Effect Size vs Number of Core Components', fontsize=14, fontweight='bold')
        
        ax.set_xlabel('Number of Core Components (k)', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add optimal value text
        if metric == 'nre':
            label_text = f'Optimal k = {best_k}\nNRE = {mean_metric_per_k[best_k]:.6f}'
        else:
            label_text = f'Optimal k = {best_k}\nMedian GMM effect size = {mean_metric_per_k[best_k]:.6f}'
        
        ax.text(0.02, 0.98, label_text, 
               transform=ax.transAxes, fontsize=11, fontweight='bold',
               bbox=dict(boxstyle="round,pad=0.4", facecolor='lightblue', alpha=0.8),
               verticalalignment='top')
        
        # Adjust layout
        plt.tight_layout()
    
    return best_k, mean_metric_per_k, all_metric_per_k, fig