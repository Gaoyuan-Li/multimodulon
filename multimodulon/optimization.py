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

if TYPE_CHECKING:
    from .core import MultiModulon

import logging
logger = logging.getLogger(__name__)


# ============================================================================
# Low-level optimization helpers
# ============================================================================

def passes_single_gene_filter(weight_vector: np.ndarray, ratio_threshold: float = 3.0) -> bool:
    """
    Determine whether a component passes the single-gene filter.
    
    Components are rejected when the largest absolute weight is more than
    ``ratio_threshold`` times the second-largest absolute weight.
    
    Args:
        weight_vector: 1D array of weights for a single component.
        ratio_threshold: Required ratio between top-1 and top-2 absolute weights.
    
    Returns:
        True if the component should be kept, False if it should be removed.
    """
    abs_weights = np.abs(weight_vector)
    
    if abs_weights.size < 2:
        # Not enough genes to assess dominance; keep the component.
        return True
    
    top_two = np.partition(abs_weights, -2)[-2:]
    largest = top_two.max()
    second_largest = top_two.min()
    
    if second_largest == 0:
        # If only one non-zero entry, treat as single-gene component (reject).
        return largest == 0

    return largest <= ratio_threshold * second_largest


def _run_single_gene_optimization(
    species_X_matrices: Dict[str, pd.DataFrame],
    k_candidates: List[int],
    num_runs: int = 1,
    mode: str = 'gpu',
    seed: int = 42,
    fig_size: Tuple[float, float] = (5, 3),
    font_path: Optional[str] = None,
    num_runs_per_dimension: Optional[int] = None
) -> Tuple[int, Dict[int, float], Dict[int, List[float]], plt.Figure]:
    """
    Optimize the number of core components using the single-gene filter metric.
    
    Components are counted when every species passes the single-gene filter
    (largest absolute weight ≥ 3 × the second-largest). Robust clustering is
    supported when ``num_runs_per_dimension`` > 1.
    """
    species_list = sorted(list(species_X_matrices.keys()))
    n_species = len(species_list)
    if n_species < 2:
        raise ValueError(f"Expected at least 2 species, got {n_species}")
    
    import torch
    import random
    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    
    if mode == 'gpu' and torch.cuda.is_available():
        print(f"Using GPU: {torch.cuda.get_device_name(0)}")
    else:
        print("Using CPU")
    
    mean_metric_per_k: Dict[int, float] = {}
    all_metric_per_k: Dict[int, List[float]] = {k: [] for k in k_candidates}
    
    for run in range(num_runs):
        k_iterator = tqdm(
            k_candidates,
            desc=f"Run {run+1}/{num_runs}",
            leave=False
        )
        
        for k in k_iterator:
            if num_runs_per_dimension is not None and num_runs_per_dimension > 1:
                from .multiview_ica import run_multiview_ica
                from sklearn.cluster import HDBSCAN
                
                a_values = {species: k for species in species_list}
                core_components = {species: [] for species in species_list}
                
                for robust_run in range(num_runs_per_dimension):
                    M_matrices = run_multiview_ica(
                        species_X_matrices,
                        a_values,
                        k,
                        mode=mode,
                        seed=seed + robust_run + run * num_runs_per_dimension
                    )
                    
                    for species in species_list:
                        M = M_matrices[species]
                        for comp_idx in range(k):
                            weight_vector = M.iloc[:, comp_idx].values
                            if passes_single_gene_filter(weight_vector):
                                component = weight_vector.copy()
                                max_idx = np.argmax(np.abs(component))
                                if component[max_idx] < 0:
                                    component *= -1
                                core_components[species].append(component)
                
                base_per_dataset = 50
                base_runs = 100
                runs_scaling_factor = num_runs_per_dimension / base_runs
                core_min_cluster_size = max(2, int(base_per_dataset * n_species * runs_scaling_factor))
                core_min_samples = max(2, int(base_per_dataset * n_species * runs_scaling_factor))
                core_min_per_dataset_count = max(1, int(base_per_dataset * runs_scaling_factor))
                
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
                    
                    clusterer = HDBSCAN(
                        min_cluster_size=core_min_cluster_size,
                        min_samples=core_min_samples,
                        cluster_selection_epsilon=0.0,
                        metric='euclidean',
                        n_jobs=-1
                    )
                    labels = clusterer.fit_predict(X)
                    
                    num_valid_components = 0
                    for label in np.unique(labels):
                        if label == -1:
                            continue
                        
                        mask = labels == label
                        cluster_size = mask.sum()
                        cluster_datasets = dataset_labels[mask]
                        dataset_counts = np.bincount(cluster_datasets, minlength=len(species_list))
                        
                        if cluster_size >= core_min_cluster_size and np.all(dataset_counts >= core_min_per_dataset_count):
                            num_valid_components += 1
                else:
                    num_valid_components = 0
            else:
                from .multiview_ica import run_multiview_ica
                
                a_values = {species: k for species in species_list}
                M_matrices = run_multiview_ica(
                    species_X_matrices,
                    a_values,
                    k,
                    mode=mode,
                    seed=seed + run
                )
                
                num_valid_components = 0
                for comp_idx in range(k):
                    if all(
                        passes_single_gene_filter(M_matrices[species].iloc[:, comp_idx].values)
                        for species in species_list
                    ):
                        num_valid_components += 1
            
            all_metric_per_k[k].append(num_valid_components)
    
    for k in k_candidates:
        if all_metric_per_k[k]:
            mean_metric_per_k[k] = np.mean(all_metric_per_k[k])
    
    k_sorted = sorted(mean_metric_per_k.keys())
    counts = [mean_metric_per_k[k] for k in k_sorted]
    
    if len(k_sorted) < 3:
        optimal_num_core_components = k_sorted[np.argmax(counts)]
    else:
        max_count = max(counts)
        min_count = min(counts)
        count_range = max_count - min_count
        
        if count_range == 0:
            optimal_num_core_components = k_sorted[0]
        else:
            threshold_percentage = 0.9
            target_count = min_count + threshold_percentage * count_range
            
            for i, (k_value, count) in enumerate(zip(k_sorted, counts)):
                if count >= target_count:
                    if i < len(k_sorted) - 2:
                        next_counts = counts[i:i+3]
                        count_std = np.std(next_counts)
                        if count_std < 0.1 * count:
                            optimal_num_core_components = k_value
                            break
                    else:
                        optimal_num_core_components = k_value
                        break
            else:
                rates = []
                for i in range(1, len(counts)):
                    rate = (counts[i] - counts[i-1]) / (k_sorted[i] - k_sorted[i-1])
                    rates.append(rate)
                
                for i in range(1, len(rates)):
                    if rates[i] < 0.1 * rates[0]:
                        optimal_num_core_components = k_sorted[i]
                        break
                else:
                    optimal_num_core_components = k_sorted[-1]

    descriptor = (
        "robust components passing filter"
        if num_runs_per_dimension and num_runs_per_dimension > 1
        else "components passing filter"
    )
    print(
        f"\nOptimal k = {optimal_num_core_components} "
        f"({descriptor} = {mean_metric_per_k[optimal_num_core_components]:.1f})"
    )

    # Create plot
    if 'plt' in globals():
        fig, ax = plt.subplots(figsize=fig_size)
        
        if font_path and os.path.exists(font_path):
            font_prop = fm.FontProperties(fname=font_path)
            plt.rcParams['font.family'] = font_prop.get_name()
    else:
        class DummyFig:
            def savefig(self, *args, **kwargs):
                pass
        fig = DummyFig()
        ax = None
    
    if ax is not None:
        k_values = sorted(mean_metric_per_k.keys())
        mean_values = [mean_metric_per_k[k] for k in k_values]
        
        ax.plot(k_values, mean_values, 'go-', linewidth=2, markersize=8,
                label='Components passing single-gene filter')
        
        ax.axvline(x=optimal_num_core_components, color='red', linestyle='--', alpha=0.7, linewidth=2,
                   label=f'Optimal k = {optimal_num_core_components}')
        ax.scatter([optimal_num_core_components], [mean_metric_per_k[optimal_num_core_components]],
                   color='red', s=120, zorder=5, edgecolors='darkred', linewidth=2)
        
        if num_runs_per_dimension is not None and num_runs_per_dimension > 1:
            ax.set_ylabel('Number of robust components', fontsize=12)
            ax.set_title('Optimization of core components (robust clustering)', fontsize=14, fontweight='bold')
        else:
            ax.set_ylabel('Number of components', fontsize=12)
            ax.set_title('Optimization of core components', fontsize=14, fontweight='bold')
        
        ax.set_xlabel('Number of Core Components (k)', fontsize=12)
        ax.legend()
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        
        if font_path and os.path.exists(font_path):
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontproperties(font_prop)  # type: ignore[name-defined]
            ax.xaxis.label.set_fontproperties(font_prop)  # type: ignore[name-defined]
            ax.yaxis.label.set_fontproperties(font_prop)  # type: ignore[name-defined]
            ax.title.set_fontproperties(font_prop)  # type: ignore[name-defined]
            for text in ax.get_legend().get_texts():
                text.set_fontproperties(font_prop)  # type: ignore[name-defined]
        
        count_desc = "robust components" if num_runs_per_dimension and num_runs_per_dimension > 1 else "components"
        label_text = (
            f'Optimal k = {optimal_num_core_components}\n'
            f'{int(round(mean_metric_per_k[optimal_num_core_components]))} {count_desc} passing filter'
        )
        
        text_obj = ax.text(
            0.02,
            0.98,
            label_text,
            transform=ax.transAxes,
            fontsize=11,
            fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.4", facecolor='lightblue', alpha=0.8),
            verticalalignment='top'
        )
        
        if font_path and os.path.exists(font_path):
            text_obj.set_fontproperties(font_prop)  # type: ignore[name-defined]
        
        plt.tight_layout()

    return optimal_num_core_components, mean_metric_per_k, all_metric_per_k, fig
# ============================================================================
# High-level optimization wrappers (from original optimization.py)
# ============================================================================

def optimize_number_of_core_components(
    multimodulon: 'MultiModulon',
    max_k: Optional[int] = None,
    step: int = 5,
    num_runs: int = 1,
    mode: str = 'gpu',
    seed: int = 42,
    save_plot: Optional[str] = None,
    save_path: Optional[str] = None,
    fig_size: Tuple[float, float] = (7, 5),
    font_path: Optional[str] = None,
    num_runs_per_dimension: Optional[int] = None
) -> int:
    """
    Optimize the number of core (shared) components across species using the
    single-gene filter metric.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        The MultiModulon instance.
    max_k : int, optional
        Maximum k to test. If None, uses the minimum dimension across species.
    step : int, default=5
        Step size for k candidates (tests k = step, 2*step, ...).
    num_runs : int, default=1
        Number of optimization runs.
    mode : str, default='gpu'
        Computation mode ('gpu' or 'cpu').
    seed : int, default=42
        Random seed for reproducibility.
    save_plot : str, optional
        Deprecated. Use ``save_path`` instead.
    save_path : str, optional
        Directory to save the optimization plot. If None, displays the plot.
    fig_size : tuple, default=(5, 3)
        Figure size as (width, height) in inches.
    font_path : str, optional
        Path to font file applied to matplotlib text elements.
    num_runs_per_dimension : int, optional
        Number of ICA runs per k value for robust clustering.
    
    Returns
    -------
    optimal_num_core_components : int
        Optimal number of core components.
    """
    species_list = list(multimodulon._species_data.keys())
    n_species = len(species_list)
    if n_species < 2:
        raise ValueError(
            f"Single-gene optimization requires at least 2 species/strains, found {n_species}"
        )
    
    print(f"Optimizing core components for {n_species} species/strains: {species_list}")
    
    for species in species_list:
        if multimodulon._species_data[species]._X is None:
            raise ValueError(
                f"X matrix not found for {species}. "
                "Please run generate_X() first to create aligned expression matrices."
            )
    
    if max_k is None:
        n_samples = [multimodulon._species_data[species].X.shape[1] for species in species_list]
        min_samples = min(n_samples)
        max_k = ((min_samples - 1) // step) * step
        max_k = max(step, max_k)
        print(f"Auto-determined max_k = {max_k} based on minimum samples ({min_samples})")
    
    k_candidates = list(range(step, max_k + 1, step))
    if not k_candidates:
        raise ValueError(f"No valid k candidates with step={step} and max_k={max_k}")
    
    species_X_matrices = {
        species: multimodulon._species_data[species].X for species in species_list
    }
    
    optimal_num_core_components, _, _, fig = _run_single_gene_optimization(
        species_X_matrices=species_X_matrices,
        k_candidates=k_candidates,
        num_runs=num_runs,
        mode=mode,
        seed=seed,
        fig_size=fig_size,
        font_path=font_path,
        num_runs_per_dimension=num_runs_per_dimension
    )
    
    if save_path or save_plot:
        if save_path:
            save_dir = Path(save_path)
            save_dir.mkdir(parents=True, exist_ok=True)
            save_file = save_dir / "num_core_optimization.svg"
        else:
            save_file = Path(save_plot)
        
        if fig is not None:
            fig.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"\nPlot saved to: {save_file}")
            plt.close(fig)
    else:
        if fig is not None:
            plt.show()
    
    multimodulon._optimal_k = optimal_num_core_components
    
    return optimal_num_core_components


def optimize_number_of_unique_components(
    multimodulon: 'MultiModulon',
    optimal_num_core_components: Optional[int] = None,
    step: int = 5,
    mode: str = 'gpu',
    seed: int = 42,
    save_plots: Optional[str] = None,
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
    save_path : str, optional
        Directory path to save the optimization plots. If None, displays plots.
        Plots are saved as 'num_unique_{species}_optimization.svg'
    fig_size : tuple, default=(5, 3)
        Figure size as (width, height) in inches
    font_path : str, optional
        Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
        If provided, uses this font for all text elements
    num_runs_per_dimension : int, default=1
        Number of ICA runs per a value for robust clustering to identify consistent components.
        Components failing the single-gene filter (largest absolute weight < 3 × second-largest)
        are discarded before clustering or counting.
        
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
                        
                        if passes_single_gene_filter(weight_vector):
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
                        
                        # Apply single-gene filter
                        weight_vector = M.iloc[:, comp_idx].values
                        
                        if passes_single_gene_filter(weight_vector):
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
