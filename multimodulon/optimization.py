"""Optimization methods for MultiModulon component analysis.

This module contains functions for optimizing the number of core and unique components
in multi-view ICA analysis, including low-level optimization functions and high-level wrappers.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Tuple, List, TYPE_CHECKING
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.font_manager as fm
import os
import time
from sklearn.model_selection import train_test_split

if TYPE_CHECKING:
    from .core import MultiModulon

import logging
logger = logging.getLogger(__name__)


# ============================================================================
# Low-level optimization functions (from multiview_ica_optimization.py)
# ============================================================================

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
    effect_size_threshold: float = 1,
    num_top_gene: int = 20,
    fig_size: Tuple[float, float] = (5, 3),
    font_path: Optional[str] = None,
    num_runs_per_dimension: Optional[int] = None
) -> Tuple[int, Dict[int, float], Dict[int, List], Dict[int, List], plt.Figure]:
    """
    Optimization using NRE or Cohen's d metric.
    
    Args:
        metric: 'nre' for Normalized Reconstruction Error, 'effect_size' for Cohen's d effect size
        num_runs_per_dimension: Number of ICA runs per k value for robust clustering (only for effect_size metric).
                               If provided, uses robust clustering to count components.
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
                if num_runs_per_dimension is not None and num_runs_per_dimension > 1:
                    # Use robust clustering approach
                    from .multiview_ica import run_multiview_ica
                    from sklearn.cluster import HDBSCAN
                    
                    # For effect_size metric, use square ICA (a=c=k)
                    a_values = {species: k for species in species_list}
                    
                    # Collect core components from multiple runs
                    core_components = {species: [] for species in species_list}
                    
                    for robust_run in range(num_runs_per_dimension):
                        # Run multi-view ICA with different seed
                        M_matrices = run_multiview_ica(
                            species_X_matrices,
                            a_values,
                            k,  # c = k (square ICA)
                            mode=mode,
                            seed=seed + robust_run + run * num_runs_per_dimension  # Different seed for each robust run
                        )
                        
                        # Process core components for each species
                        for species in species_list:
                            M = M_matrices[species]
                            for comp_idx in range(k):
                                weight_vector = M.iloc[:, comp_idx].values
                                effect_size = calculate_cohens_d_effect_size(weight_vector, seed, num_top_gene)
                                
                                if effect_size >= effect_size_threshold:
                                    # Apply sign convention
                                    component = weight_vector.copy()
                                    max_idx = np.argmax(np.abs(component))
                                    if component[max_idx] < 0:
                                        component *= -1
                                    core_components[species].append(component)
                    
                    # Cluster core components across all species
                    # Calculate scaled parameters
                    base_per_dataset = 50
                    base_runs = 100
                    runs_scaling_factor = num_runs_per_dimension / base_runs
                    core_min_cluster_size = max(2, int(base_per_dataset * n_species * runs_scaling_factor))
                    core_min_samples = max(2, int(base_per_dataset * n_species * runs_scaling_factor))
                    core_min_per_dataset_count = max(1, int(base_per_dataset * runs_scaling_factor))
                    
                    # Prepare data matrix and dataset labels
                    X_list, dataset_labels = [], []
                    for dataset_idx, species in enumerate(species_list):
                        if not core_components[species]:
                            continue
                        X = np.unique(np.array(core_components[species]), axis=0)
                        X_list.append(X)
                        dataset_labels.append(np.full(X.shape[0], dataset_idx))
                    
                    if X_list:
                        X = np.vstack(X_list)
                        dataset_labels = np.concatenate(dataset_labels)
                        
                        # Cluster components
                        clusterer = HDBSCAN(
                            min_cluster_size=core_min_cluster_size,
                            min_samples=core_min_samples,
                            cluster_selection_epsilon=0.0,
                            metric='euclidean',
                            n_jobs=-1
                        )
                        labels = clusterer.fit_predict(X)
                        
                        # Count valid clusters
                        num_above_threshold = 0
                        for label in np.unique(labels):
                            if label == -1:
                                continue
                            
                            mask = labels == label
                            cluster_size = mask.sum()
                            cluster_datasets = dataset_labels[mask]
                            dataset_counts = np.bincount(cluster_datasets, minlength=len(species_list))
                            
                            # Check if cluster meets requirements
                            if cluster_size >= core_min_cluster_size and np.all(dataset_counts >= core_min_per_dataset_count):
                                num_above_threshold += 1
                    else:
                        num_above_threshold = 0
                    
                    score = num_above_threshold
                    
                else:
                    # Original single-run approach
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
                    
                    # Count CORE components above effect_size_threshold
                    num_above_threshold = sum(1 for effect_size in core_avg_effect_sizes if effect_size > effect_size_threshold)
                    
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
        optimal_num_core_components = max(optimal_k_candidates)
    else:  # effect_size
        # More robust method: find the elbow point where we've captured most of the benefit
        k_sorted = sorted(mean_metric_per_k.keys())
        counts = [mean_metric_per_k[k] for k in k_sorted]
        
        if len(k_sorted) < 3:
            # Too few points, just take the max
            optimal_num_core_components = k_sorted[np.argmax(counts)]
        else:
            # Calculate the percentage of max count achieved at each k
            max_count = max(counts)
            min_count = min(counts)
            count_range = max_count - min_count
            
            if count_range == 0:
                # All counts are the same
                optimal_num_core_components = k_sorted[0]
            else:
                # Find where we reach 90% of the range
                threshold_percentage = 0.9
                target_count = min_count + threshold_percentage * count_range
                
                # Find the smallest k that achieves at least the target count
                for i, (k, count) in enumerate(zip(k_sorted, counts)):
                    if count >= target_count:
                        # Check if this is part of a stable plateau by looking ahead
                        if i < len(k_sorted) - 2:
                            # Check the next few points for stability
                            next_counts = counts[i:i+3]
                            count_std = np.std(next_counts)
                            # If relatively stable (std < 10% of current value), we found our k
                            if count_std < 0.1 * count:
                                optimal_num_core_components = k
                                break
                        else:
                            optimal_num_core_components = k
                            break
                else:
                    # If we didn't find a good plateau, use the elbow method
                    # Calculate the rate of change
                    if len(k_sorted) >= 3:
                        rates = []
                        for i in range(1, len(counts)):
                            rate = (counts[i] - counts[i-1]) / (k_sorted[i] - k_sorted[i-1])
                            rates.append(rate)
                        
                        # Find where the rate drops significantly (elbow)
                        for i in range(1, len(rates)):
                            if rates[i] < 0.1 * rates[0]:  # Rate dropped to less than 10% of initial rate
                                optimal_num_core_components = k_sorted[i]
                                break
                        else:
                            # If no clear elbow, take the k at 90% of max
                            optimal_num_core_components = k_sorted[-1]
                    else:
                        optimal_num_core_components = k_sorted[-1]
    
    if metric == 'nre':
        print(f"\nOptimal k = {optimal_num_core_components} (NRE = {mean_metric_per_k[optimal_num_core_components]:.6f})")
    else:
        if num_runs_per_dimension is not None and num_runs_per_dimension > 1:
            print(f"\nOptimal k = {optimal_num_core_components} (Number of robust components above threshold = {mean_metric_per_k[optimal_num_core_components]:.1f})")
        else:
            print(f"\nOptimal k = {optimal_num_core_components} (Number of components above threshold = {mean_metric_per_k[optimal_num_core_components]:.1f})")
    
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
        fig, ax = plt.subplots(figsize=fig_size)
        
        # Set font properties if provided
        if font_path and os.path.exists(font_path):
            font_prop = fm.FontProperties(fname=font_path)
            plt.rcParams['font.family'] = font_prop.get_name()
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
            ax.axvline(x=optimal_num_core_components, color='red', linestyle='--', alpha=0.7, linewidth=2,
                      label=f'Optimal k = {optimal_num_core_components}')
            ax.scatter([optimal_num_core_components], [mean_metric_per_k[optimal_num_core_components]], color='red', s=120, zorder=5,
                      edgecolors='darkred', linewidth=2)
            
            ax.set_ylabel('NRE Score (Lower is Better)', fontsize=12)
            ax.set_title('Normalized Reconstruction Error (NRE) vs Number of Core Components', fontsize=14, fontweight='bold')
            
        else:
            # Line plot for effect_size showing number of components above threshold
            mean_values = [mean_metric_per_k[k] for k in k_values]
            
            ax.plot(k_values, mean_values, 'go-', linewidth=2, markersize=8, 
                   label='QC-passed components')
            
            # Highlight optimal k
            ax.axvline(x=optimal_num_core_components, color='red', linestyle='--', alpha=0.7, linewidth=2,
                      label=f'Optimal k = {optimal_num_core_components}')
            ax.scatter([optimal_num_core_components], [mean_metric_per_k[optimal_num_core_components]], color='red', s=120, zorder=5,
                      edgecolors='darkred', linewidth=2)
            
            if num_runs_per_dimension is not None and num_runs_per_dimension > 1:
                ax.set_ylabel('Number of robust QC-passed components', fontsize=12)
                ax.set_title('Optimization of Core components (Robust clustering)', 
                            fontsize=14, fontweight='bold')
            else:
                ax.set_ylabel('Number of QC-passed components', fontsize=12)
                ax.set_title('Optimization of Core components', 
                            fontsize=14, fontweight='bold')
        
        ax.set_xlabel('Number of Core Components (k)', fontsize=12)
        ax.legend()
        
        # Apply font settings if provided
        if font_path and os.path.exists(font_path):
            # Apply to tick labels
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontproperties(font_prop)
            # Apply to axis labels
            ax.xaxis.label.set_fontproperties(font_prop)
            ax.yaxis.label.set_fontproperties(font_prop)
            # Apply to title
            ax.title.set_fontproperties(font_prop)
            # Apply to legend
            for text in ax.get_legend().get_texts():
                text.set_fontproperties(font_prop)
        
        # Force y-axis to show only integers
        if metric == 'effect_size':
            from matplotlib.ticker import MaxNLocator
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        
        # Add optimal value text
        if metric == 'nre':
            label_text = f'Optimal k = {optimal_num_core_components}\nNRE = {mean_metric_per_k[optimal_num_core_components]:.6f}'
        else:
            if num_runs_per_dimension is not None and num_runs_per_dimension > 1:
                label_text = f'Optimal k = {optimal_num_core_components}\n{int(mean_metric_per_k[optimal_num_core_components])} robust components above threshold'
            else:
                label_text = f'Optimal k = {optimal_num_core_components}\n{int(mean_metric_per_k[optimal_num_core_components])} components above threshold'
        
        text_obj = ax.text(0.02, 0.98, label_text, 
               transform=ax.transAxes, fontsize=11, fontweight='bold',
               bbox=dict(boxstyle="round,pad=0.4", facecolor='lightblue', alpha=0.8),
               verticalalignment='top')
        
        # Apply font to text box if provided
        if font_path and os.path.exists(font_path):
            text_obj.set_fontproperties(font_prop)
        
        # Adjust layout
        plt.tight_layout()
    
    return optimal_num_core_components, mean_metric_per_k, all_metric_per_k, fig


# ============================================================================
# High-level optimization wrappers (from original optimization.py)
# ============================================================================

def optimize_number_of_core_components(
    multimodulon: 'MultiModulon',
    max_k: Optional[int] = None,
    step: int = 5,
    max_a_per_view: Optional[int] = None,
    train_frac: float = 0.75,
    num_runs: int = 1,
    mode: str = 'gpu',
    seed: int = 42,
    save_plot: Optional[str] = None,
    metric: str = 'effect_size',
    threshold: Optional[float] = None,
    effect_size_threshold: float = 1,
    num_top_gene: int = 20,
    save_path: Optional[str] = None,
    fig_size: Tuple[float, float] = (7, 5),
    font_path: Optional[str] = None,
    num_runs_per_dimension: Optional[int] = None
) -> Tuple[int, Dict[int, float]]:
    """
    Optimize the number of core components using specified metric.
    
    Available Metrics
    -----------------
    'nre' : Normalized Reconstruction Error
        Measures component sharing across views. Minimum at optimal k.
        
    'effect_size' : Cohen's d Effect Size  
        Measures gene weight separation (top genes vs rest).
        Higher values indicate clearer component structure.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        The MultiModulon instance
    max_k : int, optional
        Maximum k to test. If None, uses the maximum dimension from generate_X
    step : int, default=5
        Step size for k candidates (tests k = step, 2*step, 3*step, ...)
    max_a_per_view : int, optional
        Maximum components per view. If None, uses max_k
    train_frac : float, default=0.75
        Fraction of data to use for training
    num_runs : int, default=1
        Number of runs for optimization. For 'effect_size' metric, this controls
        the number of ICA runs used for robust clustering to identify components
    mode : str, default='gpu'
        'gpu' or 'cpu' mode
    seed : int, default=42
        Random seed for reproducibility
    save_plot : str, optional
        Path to save the metric vs k plot. If None, displays the plot
    metric : str, default='nre'
        Optimization metric: 'nre' or 'effect_size' (Cohen's d effect size)
    threshold : float, optional
        For effect_size metric: Cohen's d threshold for effect size analysis and visualization
    effect_size_threshold : float, default=1
        Minimum Cohen's d effect size threshold for counting components.
        Only components with effect size above this threshold are counted as "effective".
    num_top_gene : int, default=20
        Number of top genes to use when calculating Cohen's d effect size
    save_path : str, optional
        Directory path to save the optimization plot. If None, displays plot.
        Plot is saved as 'num_core_optimization.svg'
    fig_size : tuple, default=(5, 3)
        Figure size as (width, height) in inches
    font_path : str, optional
        Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
        If provided, uses this font for all text elements
    num_runs_per_dimension : int, optional
        Number of ICA runs per k value for robust clustering (only for 'effect_size' metric).
        If provided, uses robust clustering to identify consistent components
    
    Returns
    -------
    For metric='effect_size':
        optimal_num_core_components : int
            Optimal number of core components
    
    For metric='nre':
        optimal_num_core_components : int
            Optimal number of core components
        metric_scores : dict
            Dictionary mapping k values to mean NRE scores
    
    Examples
    --------
    >>> # Using effect size metric
    >>> optimal_num_core_components, scores = multiModulon.optimize_number_of_core_components(
    ...     metric='effect_size',
    ...     effect_size_threshold=5
    ... )
    """
    # Check prerequisites
    species_list = list(multimodulon._species_data.keys())
    n_species = len(species_list)
    if n_species < 2:
        raise ValueError(
            f"NRE optimization requires at least 2 species/strains, found {n_species}"
        )
    
    print(f"Optimizing core components for {n_species} species/strains: {species_list}")
    
    # Check if all species have X matrices
    for species in species_list:
        if multimodulon._species_data[species]._X is None:
            raise ValueError(
                f"X matrix not found for {species}. "
                "Please run generate_X() first to create aligned expression matrices."
            )
    
    # Determine max_k if not provided
    if max_k is None:
        # Get number of samples for each species (should be consistent after alignment)
        n_samples = []
        for species in species_list:
            X = multimodulon._species_data[species].X
            n_samples.append(X.shape[1])  # number of samples (columns)
        min_samples = min(n_samples)
        # Set max_k to the largest multiple of step that's less than min_samples
        # This ensures k candidates don't exceed data constraints
        max_k = ((min_samples - 1) // step) * step
        max_k = max(step, max_k)  # Ensure at least one k candidate
        print(f"Auto-determined max_k = {max_k} based on minimum samples ({min_samples})")
    
    # Generate k candidates
    k_candidates = list(range(step, max_k + 1, step))
    if not k_candidates:
        raise ValueError(f"No valid k candidates with step={step} and max_k={max_k}")
    
    # Set max_a_per_view if not provided
    if max_a_per_view is None:
        max_a_per_view = max_k
    
    # Prepare X matrices
    species_X_matrices = {}
    for species in species_list:
        species_X_matrices[species] = multimodulon._species_data[species].X
    
    # Run optimization
    optimal_num_core_components, metric_scores, all_metric_per_k, fig = run_nre_optimization(
        species_X_matrices=species_X_matrices,
        k_candidates=k_candidates,
        max_a_per_view=max_a_per_view,
        train_frac=train_frac,
        num_runs=num_runs,
        mode=mode,
        seed=seed,
        metric=metric,
        threshold=threshold,
        effect_size_threshold=effect_size_threshold,
        num_top_gene=num_top_gene,
        fig_size=fig_size,
        font_path=font_path,
        num_runs_per_dimension=num_runs_per_dimension
    )
    
    # Handle save_plot (deprecated) and save_path
    if save_path or save_plot:
        # Determine save file path
        if save_path:
            save_path = Path(save_path)
            save_path.mkdir(parents=True, exist_ok=True)
            save_file = save_path / "num_core_optimization.svg"
        else:
            # Use deprecated save_plot parameter
            save_file = Path(save_plot)
        
        fig.savefig(save_file, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved to: {save_file}")
        plt.close(fig)
    else:
        plt.show()
    
    # Store the optimal k for later use
    multimodulon._optimal_k = optimal_num_core_components
    
    # Return different values based on metric
    if metric == 'effect_size':
        return optimal_num_core_components
    else:  # metric == 'nre'
        return optimal_num_core_components, metric_scores


def optimize_number_of_unique_components(
    multimodulon: 'MultiModulon',
    optimal_num_core_components: Optional[int] = None,
    step: int = 5,
    mode: str = 'gpu',
    seed: int = 42,
    save_plots: Optional[str] = None,
    effect_size_threshold: float = 1,
    num_top_gene: int = 20,
    save_path: Optional[str] = None,
    fig_size: Tuple[float, float] = (7, 5),
    font_path: Optional[str] = None,
    num_runs_per_dimension: int = 1
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Optimize the number of unique components for each species.
    
    This method runs ICA with varying numbers of unique components for each species
    while keeping other species at c=optimal_num_core_components, and finds the optimal number that
    maximizes consistent unique components.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        The MultiModulon instance
    optimal_num_core_components : int, optional
        Number of core components. If None, uses the value from optimize_number_of_core_components
    step : int, default=5
        Step size for testing unique components (tests optimal_num_core_components+5, optimal_num_core_components+10, ...)
    mode : str, default='gpu'
        'gpu' or 'cpu' mode
    seed : int, default=42
        Random seed for reproducibility
    save_plots : str, optional
        Directory to save plots. If None, displays the plots
    effect_size_threshold : float, default=1
        Minimum Cohen's d effect size threshold for counting components
    num_top_gene : int, default=20
        Number of top genes to use when calculating Cohen's d effect size
    save_path : str, optional
        Directory path to save the optimization plots. If None, displays plots.
        Plots are saved as 'num_unique_{species}_optimization.svg'
    fig_size : tuple, default=(5, 3)
        Figure size as (width, height) in inches
    font_path : str, optional
        Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
        If provided, uses this font for all text elements
    num_runs_per_dimension : int, default=1
        Number of ICA runs per a value for robust clustering to identify consistent components
        
    Returns
    -------
    optimal_num_unique_components : dict
        Dictionary mapping species names to optimal number of unique components
    optimal_num_total_components : dict
        Dictionary mapping species names to optimal total components (a values).
        Can be used directly with run_robust_multiview_ica.
    """
    # Use stored optimal k if not provided
    if optimal_num_core_components is None:
        if hasattr(multimodulon, '_optimal_k'):
            optimal_num_core_components = multimodulon._optimal_k
        else:
            raise ValueError("optimal_num_core_components not provided and no optimal k found. Run optimize_number_of_core_components first.")
    
    print(f"\nOptimizing unique components with core k = {optimal_num_core_components}")
    
    species_list = list(multimodulon._species_data.keys())
    n_species = len(species_list)
    
    # Prepare X matrices
    species_X_matrices = {}
    for species in species_list:
        species_X_matrices[species] = multimodulon._species_data[species].X
    
    # Results storage
    optimal_num_unique_components = {}
    
    # Process each species
    for target_species in species_list:
        print(f"\n{'='*60}")
        print(f"Optimizing unique components for {target_species}")
        print(f"{'='*60}")
        
        # Get minimum dimension (number of samples)
        min_dim = species_X_matrices[target_species].shape[1]
        
        # Generate a candidates: from optimal_num_core_components to closest multiple of 5 below min_dim
        a_candidates = []
        a = optimal_num_core_components
        while a < min_dim - 5:
            a_candidates.append(a)
            a += step
        # Add the last valid a
        if a_candidates and a_candidates[-1] < min_dim - 5:
            last_valid = ((min_dim - 5) // step) * step
            if last_valid not in a_candidates:
                a_candidates.append(last_valid)
        
        if not a_candidates:
            a_candidates = [optimal_num_core_components]
        
        print(f"Testing a values: {a_candidates}")
        
        # Store results
        consistent_counts = {}
        
        # Test each a value
        for a_test in tqdm(a_candidates, desc=f"Testing a values for {target_species}"):
            if num_runs_per_dimension > 1:
                # Use robust clustering approach
                from sklearn.cluster import HDBSCAN
                
                # Set up a_values: a_test for target species, optimal_num_core_components for others
                a_values = {}
                for species in species_list:
                    if species == target_species:
                        a_values[species] = a_test
                    else:
                        a_values[species] = optimal_num_core_components
                
                # Collect unique components from multiple runs
                unique_components_runs = []
                
                for run_idx in range(num_runs_per_dimension):
                    # Run ICA
                    from .multiview_ica import run_multiview_ica
                    M_matrices = run_multiview_ica(
                        species_X_matrices,
                        a_values,
                        optimal_num_core_components,  # c = optimal_num_core_components
                        mode=mode,
                        seed=seed + run_idx  # Different seed for each run
                    )
                    
                    # Extract unique components for target species
                    M = M_matrices[target_species]
                    for comp_idx in range(optimal_num_core_components, a_test):
                        weight_vector = M.iloc[:, comp_idx].values
                        effect_size = calculate_cohens_d_effect_size(weight_vector, seed, num_top_gene)
                        
                        if effect_size >= effect_size_threshold:
                            # Apply sign convention
                            component = weight_vector.copy()
                            max_idx = np.argmax(np.abs(component))
                            if component[max_idx] < 0:
                                component *= -1
                            unique_components_runs.append(component)
                
                # Cluster unique components
                if unique_components_runs:
                    X = np.unique(np.array(unique_components_runs), axis=0)
                    
                    # Calculate scaled parameters
                    base_unique_min = 50
                    base_runs = 100
                    runs_scaling_factor = num_runs_per_dimension / base_runs
                    unique_min_cluster_size = max(2, int(base_unique_min * runs_scaling_factor))
                    unique_min_samples = max(2, int(base_unique_min * runs_scaling_factor))
                    
                    # Cluster components
                    clusterer = HDBSCAN(
                        min_cluster_size=unique_min_cluster_size,
                        min_samples=unique_min_samples,
                        cluster_selection_epsilon=0.0,
                        metric='euclidean',
                        n_jobs=-1
                    )
                    labels = clusterer.fit_predict(X)
                    
                    # Count valid clusters
                    consistent_unique = 0
                    for label in np.unique(labels):
                        if label == -1:
                            continue
                        mask = labels == label
                        cluster_size = mask.sum()
                        if cluster_size >= unique_min_cluster_size:
                            consistent_unique += 1
                else:
                    consistent_unique = 0
                
                consistent_counts[a_test] = consistent_unique
                
            else:
                # Original single-run approach
                # Set up a_values: a_test for target species, optimal_num_core_components for others
                a_values = {}
                for species in species_list:
                    if species == target_species:
                        a_values[species] = a_test
                    else:
                        a_values[species] = optimal_num_core_components
                
                # Run ICA
                from .multiview_ica import run_multiview_ica
                M_matrices = run_multiview_ica(
                    species_X_matrices,
                    a_values,
                    optimal_num_core_components,  # c = optimal_num_core_components
                    mode=mode,
                    seed=seed
                )
                
                # Generate A matrix for target species
                M = M_matrices[target_species]
                X = species_X_matrices[target_species]
                A = M.T @ X
                A.index = M.columns
                A.columns = X.columns
                
                # Get sample sheet
                sample_sheet = multimodulon._species_data[target_species].sample_sheet
                
                # Count consistent unique components
                consistent_unique = 0
                
                for comp in A.index:
                    if comp.startswith("Unique"):
                        # Get component index (Unique components start after core components)
                        comp_idx = int(comp.split('_')[1]) - 1 + optimal_num_core_components
                        
                        # Check Cohen's d threshold
                        weight_vector = M.iloc[:, comp_idx].values
                        effect_size = calculate_cohens_d_effect_size(weight_vector, seed, num_top_gene)
                        
                        if effect_size >= effect_size_threshold:
                            # Check consistency
                            if _check_component_consistency(multimodulon, A, sample_sheet, comp):
                                consistent_unique += 1
                
                consistent_counts[a_test] = consistent_unique
        
        # Find optimal a using 90% growth method
        a_sorted = sorted(consistent_counts.keys())
        counts = [consistent_counts[a] for a in a_sorted]
        
        if len(a_sorted) < 2:
            optimal_a = a_sorted[0]
        else:
            max_count = max(counts)
            min_count = min(counts)
            count_range = max_count - min_count
            
            if count_range == 0:
                optimal_a = a_sorted[0]
            else:
                # Find where we reach 90% of the range
                threshold_percentage = 0.9
                target_count = min_count + threshold_percentage * count_range
                
                # Find the smallest a that achieves at least the target count
                optimal_a = a_sorted[0]
                for i, (a, count) in enumerate(zip(a_sorted, counts)):
                    if count >= target_count:
                        optimal_a = a
                        break
        
        optimal_num_unique_components[target_species] = optimal_a - optimal_num_core_components
        
        # Create plot
        fig, ax = plt.subplots(figsize=fig_size)
        
        # Set font properties if provided
        if font_path and os.path.exists(font_path):
            font_prop = fm.FontProperties(fname=font_path)
            plt.rcParams['font.family'] = font_prop.get_name()
        
        ax.plot(a_sorted, counts, 'bo-', linewidth=2, markersize=8, 
               label='QC-passed components')
        
        # Highlight optimal a
        ax.axvline(x=optimal_a, color='red', linestyle='--', alpha=0.7, linewidth=2,
                  label=f'Optimal a = {optimal_a}')
        ax.scatter([optimal_a], [consistent_counts[optimal_a]], color='red', s=120, zorder=5,
                  edgecolors='darkred', linewidth=2)
        
        ax.set_xlabel('Number of Total Components (a)', fontsize=12)
        if num_runs_per_dimension > 1:
            ax.set_ylabel('Number of robust QC-passed components', fontsize=12)
            ax.set_title(f'Optimization of Unique Components for {target_species} (Robust clustering)', fontsize=14, fontweight='bold')
        else:
            ax.set_ylabel('Number of QC-passed components', fontsize=12)
            ax.set_title(f'Optimization of Unique Components for {target_species}', fontsize=14, fontweight='bold')
        ax.legend()
        
        # Set font for tick labels if font_path provided
        if font_path and os.path.exists(font_path):
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontproperties(font_prop)
            ax.xaxis.label.set_fontproperties(font_prop)
            ax.yaxis.label.set_fontproperties(font_prop)
            ax.title.set_fontproperties(font_prop)
            # Update legend font
            for text in ax.get_legend().get_texts():
                text.set_fontproperties(font_prop)
        
        # Force integer y-axis
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        
        # Add text box with result
        if num_runs_per_dimension > 1:
            label_text = f'Optimal a = {optimal_a}\n{consistent_counts[optimal_a]} robust unique components'
        else:
            label_text = f'Optimal a = {optimal_a}\n{consistent_counts[optimal_a]} consistent unique components'
        text_obj = ax.text(0.02, 0.98, label_text, 
               transform=ax.transAxes, fontsize=11, fontweight='bold',
               bbox=dict(boxstyle="round,pad=0.4", facecolor='lightblue', alpha=0.8),
               verticalalignment='top')
        
        # Apply font to text box if font_path provided
        if font_path and os.path.exists(font_path):
            text_obj.set_fontproperties(font_prop)
        
        plt.tight_layout()
        
        # Handle save_plots (deprecated) and save_path
        if save_path or save_plots:
            # Determine save directory
            if save_path:
                save_dir = Path(save_path)
            else:
                # Use deprecated save_plots parameter
                save_dir = Path(save_plots)
            
            save_dir.mkdir(parents=True, exist_ok=True)
            plot_file = save_dir / f"num_unique_{target_species}_optimization.svg"
            fig.savefig(plot_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to: {plot_file}")
            plt.close(fig)
        else:
            plt.show()
        
        if num_runs_per_dimension > 1:
            print(f"\nOptimal a for {target_species}: {optimal_a} ({consistent_counts[optimal_a]} robust unique components)")
        else:
            print(f"\nOptimal a for {target_species}: {optimal_a} ({consistent_counts[optimal_a]} consistent unique components)")
    
    # Create optimal total components dictionary
    optimal_num_total_components = {}
    for species, num_unique in optimal_num_unique_components.items():
        optimal_num_total_components[species] = num_unique + optimal_num_core_components
    
    print(f"\n{'='*60}")
    print("Optimization Summary")
    print(f"{'='*60}")
    print(f"Core components (c): {optimal_num_core_components}")
    for species, num_unique in optimal_num_unique_components.items():
        print(f"{species}: a = {optimal_num_total_components[species]} (unique components: {num_unique})")
    
    return optimal_num_unique_components, optimal_num_total_components


# ============================================================================
# Helper functions
# ============================================================================

def _check_component_consistency(multimodulon: 'MultiModulon', A: pd.DataFrame, sample_sheet: pd.DataFrame, component_name: str) -> bool:
    """
    Check if a component is consistent across biological replicates.
    
    Returns True if component is consistent, False if it should be dropped.
    """
    # Build replicate groups
    if "biological_replicate" not in sample_sheet.columns:
        # If no biological_replicate column, assume all samples are independent
        return True
    
    st = sample_sheet.copy()
    grp_id = 0
    replicate_groups = []
    for v in st["biological_replicate"]:
        if v == 1 and replicate_groups:
            grp_id += 1
        replicate_groups.append(grp_id)
    st["replicate_group"] = replicate_groups
    sample_to_grp = st["replicate_group"].to_dict()
    
    # Get component values
    row = A.loc[component_name]
    
    # Check each replicate group
    for gid in st["replicate_group"].unique():
        cols = [c for c in A.columns if c in sample_to_grp and sample_to_grp[c] == gid]
        if not cols:
            continue
            
        vals = row[cols].astype(float).values
        if vals.size < 2:
            continue
        
        # Check pairwise consistency between all biological replicates
        for i in range(len(vals)):
            for j in range(i + 1, len(vals)):
                val1, val2 = vals[i], vals[j]
                abs1, abs2 = abs(val1), abs(val2)
                
                # Check if this pair has large opposite signs
                opposite_signs = (val1 > 0 and val2 < 0) or (val1 < 0 and val2 > 0)
                both_large = abs1 > 5 and abs2 > 5
                
                if opposite_signs and both_large:
                    return False  # Inconsistent
                
                # Check if this pair has wide span
                pair_span = abs(abs1 - abs2) > 20
                if pair_span:
                    return False  # Inconsistent
    
    return True  # Consistent