"""Optimization methods for MultiModulon component analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Tuple, TYPE_CHECKING
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from .multiview_ica import run_multiview_ica
from .multiview_ica_optimization import run_nre_optimization, calculate_cohens_d_effect_size

if TYPE_CHECKING:
    from .core import MultiModulon

import logging
logger = logging.getLogger(__name__)


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
    metric: str = 'nre',
    threshold: Optional[float] = None,
    effect_size_threshold: float = 5,
    num_top_gene: int = 20
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
        Number of cross-validation runs
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
    effect_size_threshold : float, default=5
        Minimum Cohen's d effect size threshold for counting components.
        Only components with effect size above this threshold are counted as "effective".
    num_top_gene : int, default=20
        Number of top genes to use when calculating Cohen's d effect size
    
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
    k_candidates = list(range(step, min(max_k + 1, 100), step))
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
        num_top_gene=num_top_gene
    )
    
    # Save or display plot
    if save_plot:
        fig.savefig(save_plot, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved to: {save_plot}")
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
    effect_size_threshold: float = 5,
    num_top_gene: int = 20
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
    effect_size_threshold : float, default=5
        Minimum Cohen's d effect size threshold for counting components
    num_top_gene : int, default=20
        Number of top genes to use when calculating Cohen's d effect size
        
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
            # Set up a_values: a_test for target species, optimal_num_core_components for others
            a_values = {}
            for species in species_list:
                if species == target_species:
                    a_values[species] = a_test
                else:
                    a_values[species] = optimal_num_core_components
            
            # Run ICA
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
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(a_sorted, counts, 'bo-', linewidth=2, markersize=8, 
               label='Consistent unique components')
        
        # Highlight optimal a
        ax.axvline(x=optimal_a, color='red', linestyle='--', alpha=0.7, linewidth=2,
                  label=f'Optimal a = {optimal_a}')
        ax.scatter([optimal_a], [consistent_counts[optimal_a]], color='red', s=120, zorder=5,
                  edgecolors='darkred', linewidth=2)
        
        ax.set_xlabel('Number of Total Components (a)', fontsize=12)
        ax.set_ylabel('Number of Consistent Unique Components', fontsize=12)
        ax.set_title(f'Optimization of Unique Components for {target_species}', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Force integer y-axis
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        
        # Add text box with result
        label_text = f'Optimal a = {optimal_a}\n{consistent_counts[optimal_a]} consistent unique components'
        ax.text(0.02, 0.98, label_text, 
               transform=ax.transAxes, fontsize=11, fontweight='bold',
               bbox=dict(boxstyle="round,pad=0.4", facecolor='lightblue', alpha=0.8),
               verticalalignment='top')
        
        plt.tight_layout()
        
        # Save or display plot
        if save_plots:
            Path(save_plots).mkdir(exist_ok=True)
            plot_path = Path(save_plots) / f"unique_optimization_{target_species}.png"
            fig.savefig(plot_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to: {plot_path}")
        else:
            plt.show()
        
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
        if vals.size == 0:
            continue
        
        abs_vals = np.abs(vals)
        mixed_signs = (vals > 0).any() and (vals < 0).any()
        all_big = (abs_vals > 5).all()
        wide_span = abs_vals.max() - abs_vals.min() > 20
        
        if (mixed_signs and all_big) or wide_span:
            return False  # Inconsistent
    
    return True  # Consistent