"""Core iModulon stability analysis functions."""

from __future__ import annotations
from typing import Dict, Optional, Tuple, List, Union
import pandas as pd
import numpy as np
import logging
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch
from scipy.stats import pearsonr, median_abs_deviation
from scipy.signal import argrelextrema
from sklearn.cluster import DBSCAN

logger = logging.getLogger(__name__)


def core_iModulon_stability(
    multimodulon,
    component: str,
    threshold_method: str = "mad",
    manual_threshold: Optional[float] = None,
    save_path: Optional[str] = None,
    fig_size: Tuple[float, float] = (10, 6),
    font_path: Optional[str] = None,
    show_stats: bool = True
) -> Tuple[List[str], float, Dict[str, float]]:
    """
    Quantify core iModulon stability across species using pairwise correlations.
    
    This function calculates how similar a core iModulon is across different species
    by computing the mean pairwise correlation of M matrix weights. It automatically
    determines a threshold to classify species as stable or unstable.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        MultiModulon instance containing the data
    component : str
        Core component name (e.g., 'Core_1')
    threshold_method : str, optional
        Method for determining stability threshold:
        - "mad": Modified Z-score with Median Absolute Deviation (default)
        - "clustering": DBSCAN clustering to identify groups
        - "otsu": Otsu's method for optimal threshold
        - "elbow": Find elbow point in sorted scores
        - "manual": Use user-specified threshold
    manual_threshold : float, optional
        Threshold value when threshold_method="manual"
    save_path : str, optional
        Path to save the plot. If None, displays the plot
    fig_size : tuple, optional
        Figure size as (width, height). Default: (10, 6)
    font_path : str, optional
        Path to custom font file
    show_stats : bool, optional
        Whether to show statistics on the plot. Default: True
        
    Returns
    -------
    stable_species : list
        List of species names classified as stable
    threshold : float
        The threshold value used for classification
    stability_scores : dict
        Dictionary mapping species names to stability scores
        
    Raises
    ------
    ValueError
        If component not found, insufficient species, or invalid parameters
        
    Examples
    --------
    >>> # Automatic threshold detection using MAD
    >>> stable, thresh, scores = mm.core_iModulon_stability("Core_1")
    
    >>> # Manual threshold specification
    >>> stable, thresh, scores = mm.core_iModulon_stability(
    ...     "Core_1", 
    ...     threshold_method="manual", 
    ...     manual_threshold=0.7
    ... )
    """
    # Validate inputs
    
    # Find species that have this component
    species_with_component = []
    component_weights = {}
    
    for species_name in multimodulon.species:
        species_data = multimodulon[species_name]
        if species_data._M is not None and component in species_data.M.columns:
            species_with_component.append(species_name)
            # Extract weights for this component
            component_weights[species_name] = species_data.M[component].values
    
    if len(species_with_component) < 2:
        raise ValueError(f"Component '{component}' found in fewer than 2 species. "
                        f"Found in: {species_with_component}")
    
    # Calculate pairwise correlations
    n_species = len(species_with_component)
    correlation_matrix = np.ones((n_species, n_species))
    
    for i, species1 in enumerate(species_with_component):
        for j, species2 in enumerate(species_with_component):
            if i < j:
                corr, _ = pearsonr(component_weights[species1], component_weights[species2])
                correlation_matrix[i, j] = corr
                correlation_matrix[j, i] = corr
    
    # Calculate mean correlation for each species (excluding self-correlation)
    stability_scores = {}
    for i, species in enumerate(species_with_component):
        # Mean of correlations with other species
        other_correlations = np.concatenate([correlation_matrix[i, :i], correlation_matrix[i, i+1:]])
        stability_scores[species] = np.mean(other_correlations)
    
    # Determine threshold
    scores_array = np.array(list(stability_scores.values()))
    
    if threshold_method == "manual":
        if manual_threshold is None:
            raise ValueError("manual_threshold must be provided when threshold_method='manual'")
        threshold = manual_threshold
        
    elif threshold_method == "mad":
        threshold = _mad_threshold(scores_array)
        
    elif threshold_method == "clustering":
        threshold = _clustering_threshold(scores_array)
        
    elif threshold_method == "otsu":
        threshold = _otsu_threshold(scores_array)
        
    elif threshold_method == "elbow":
        threshold = _elbow_threshold(scores_array)
        
    else:
        raise ValueError(f"Unknown threshold_method: {threshold_method}")
    
    # Classify species
    stable_species = [species for species, score in stability_scores.items() 
                     if score >= threshold]
    
    # Create visualization
    _plot_stability(stability_scores, threshold, component, species_with_component,
                   stable_species, save_path, fig_size, font_path, show_stats,
                   threshold_method)
    
    return stable_species, threshold, stability_scores


def _mad_threshold(scores: np.ndarray) -> float:
    """
    Calculate threshold using Modified Z-score with Median Absolute Deviation.
    
    Only considers species unstable when they are truly significant outliers.
    Uses a more conservative approach to avoid unnecessary unstable classifications.
    """
    if len(scores) < 3:
        # Too few points, consider all stable
        return np.min(scores) - 0.01
    
    median = np.median(scores)
    mad = median_abs_deviation(scores)
    score_range = np.max(scores) - np.min(scores)
    
    # If MAD is very small or range is small, consider all species stable
    if mad < 0.05 or score_range < 0.1:
        return np.min(scores) - 0.01  # All species are stable
    
    # Modified Z-scores with more stringent threshold
    modified_z_scores = 0.6745 * (scores - median) / mad
    
    # Use stricter threshold (-3.0 instead of -2.5) to be more conservative
    # Only mark as outlier if it's a very clear deviation
    outlier_mask = modified_z_scores < -3.0
    
    if np.any(outlier_mask):
        # Only create threshold if outliers are significantly different
        non_outliers = scores[~outlier_mask]
        outliers = scores[outlier_mask]
        gap = np.min(non_outliers) - np.max(outliers)
        
        # If gap between outliers and non-outliers is small, consider all stable
        if gap < 0.1:
            return np.min(scores) - 0.01
        
        threshold = (np.max(outliers) + np.min(non_outliers)) / 2
    else:
        # No clear outliers, all species are stable
        return np.min(scores) - 0.01
    
    return threshold


def _clustering_threshold(scores: np.ndarray) -> float:
    """
    Use DBSCAN clustering to identify groups and set threshold between them.
    Only creates unstable groups when there's a significant separation.
    """
    if len(scores) < 3:
        # Too few points, consider all stable
        return np.min(scores) - 0.01
    
    # Reshape for sklearn
    X = scores.reshape(-1, 1)
    
    # Adaptive eps based on score range
    score_range = np.max(scores) - np.min(scores)
    if score_range < 0.15:  # More lenient threshold
        # Very similar scores, consider all stable
        return np.min(scores) - 0.01
    
    # Try different eps values
    best_eps = None
    best_n_clusters = 0
    
    for eps in np.linspace(score_range * 0.1, score_range * 0.5, 10):
        clusterer = DBSCAN(eps=eps, min_samples=1)
        labels = clusterer.fit_predict(X)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        
        # Prefer 2 clusters
        if n_clusters == 2:
            best_eps = eps
            best_n_clusters = 2
            break
        elif n_clusters > 1 and n_clusters < len(scores):
            if best_n_clusters == 0 or abs(n_clusters - 2) < abs(best_n_clusters - 2):
                best_eps = eps
                best_n_clusters = n_clusters
    
    if best_eps is None:
        # Clustering failed, fall back to MAD
        return _mad_threshold(scores)
    
    # Get final clustering
    clusterer = DBSCAN(eps=best_eps, min_samples=1)
    labels = clusterer.fit_predict(X)
    
    # Find threshold between clusters
    unique_labels = np.unique(labels[labels != -1])
    if len(unique_labels) < 2:
        return _mad_threshold(scores)
    
    # Get cluster means and check separation
    cluster_means = []
    cluster_sizes = []
    for label in unique_labels:
        cluster_scores = scores[labels == label]
        cluster_means.append(np.mean(cluster_scores))
        cluster_sizes.append(len(cluster_scores))
    
    cluster_means = np.array(cluster_means)
    cluster_sizes = np.array(cluster_sizes)
    
    # Sort by cluster means
    sorted_indices = np.argsort(cluster_means)
    cluster_means = cluster_means[sorted_indices]
    cluster_sizes = cluster_sizes[sorted_indices]
    
    # Check if there's significant separation between clusters
    if len(cluster_means) >= 2:
        gap = cluster_means[1] - cluster_means[0]
        # Only create threshold if gap is substantial
        if gap < 0.1:
            return np.min(scores) - 0.01  # All stable
        
        threshold = (cluster_means[0] + cluster_means[1]) / 2
    else:
        # Single cluster, all stable
        return np.min(scores) - 0.01
    
    return threshold


def _otsu_threshold(scores: np.ndarray) -> float:
    """
    Apply Otsu's method to find optimal threshold.
    Only creates threshold if variance improvement is substantial.
    """
    if len(scores) < 3:
        # Too few points, consider all stable
        return np.min(scores) - 0.01
    
    # If scores are very similar, consider all stable
    if np.std(scores) < 0.05 or (np.max(scores) - np.min(scores)) < 0.1:
        return np.min(scores) - 0.01
    
    # Create histogram
    n_bins = min(10, len(scores))
    counts, bin_edges = np.histogram(scores, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Otsu algorithm
    best_threshold = np.mean(scores)
    max_variance = 0
    
    total = counts.sum()
    if total == 0:
        return np.mean(scores)
    
    sum_total = (counts * bin_centers).sum()
    
    w0 = 0
    sum0 = 0
    
    for i in range(len(counts) - 1):
        w0 += counts[i]
        if w0 == 0:
            continue
        
        w1 = total - w0
        if w1 == 0:
            break
        
        sum0 += counts[i] * bin_centers[i]
        m0 = sum0 / w0
        m1 = (sum_total - sum0) / w1
        
        # Between-class variance
        variance = w0 * w1 * (m0 - m1) ** 2
        
        if variance > max_variance:
            max_variance = variance
            best_threshold = bin_centers[i]
    
    # Only use threshold if the variance improvement is meaningful
    # Check if the separation is significant enough
    total_variance = np.var(scores) * len(scores)
    if max_variance < total_variance * 0.1:  # Less than 10% improvement
        return np.min(scores) - 0.01  # All stable
    
    # Also check if threshold creates meaningful separation
    stable_scores = scores[scores >= best_threshold]
    unstable_scores = scores[scores < best_threshold]
    
    if len(unstable_scores) == 0:
        return np.min(scores) - 0.01  # All stable
    
    # Check gap between groups
    if len(stable_scores) > 0:
        gap = np.min(stable_scores) - np.max(unstable_scores)
        if gap < 0.05:  # Very small gap
            return np.min(scores) - 0.01  # All stable
    
    return best_threshold


def _elbow_threshold(scores: np.ndarray) -> float:
    """
    Find elbow point in sorted scores for threshold.
    Only creates threshold if elbow represents significant change.
    """
    if len(scores) < 3:
        # Too few points, consider all stable
        return np.min(scores) - 0.01
    
    # Sort scores
    sorted_scores = np.sort(scores)
    
    # If scores are very similar or range is small
    if np.std(sorted_scores) < 0.05 or (np.max(scores) - np.min(scores)) < 0.1:
        return np.min(sorted_scores) - 0.01
    
    # Create line from first to last point
    n_points = len(sorted_scores)
    x = np.arange(n_points)
    
    # Calculate perpendicular distance from each point to the line
    line_vec = np.array([n_points - 1, sorted_scores[-1] - sorted_scores[0]])
    line_vec = line_vec / np.linalg.norm(line_vec)
    
    distances = []
    for i in range(n_points):
        point_vec = np.array([i, sorted_scores[i] - sorted_scores[0]])
        # Perpendicular distance
        dist = np.abs(np.cross(line_vec, point_vec))
        distances.append(dist)
    
    # Find maximum distance (elbow point)
    max_distance = np.max(distances)
    elbow_idx = np.argmax(distances)
    
    # Check if elbow is significant enough
    score_range = sorted_scores[-1] - sorted_scores[0]
    if max_distance < score_range * 0.1:  # Less than 10% of range
        return np.min(scores) - 0.01  # All stable
    
    if elbow_idx == 0 or elbow_idx == n_points - 1:
        # Elbow at extremes, consider all stable
        return np.min(scores) - 0.01
    else:
        # Check if elbow creates meaningful separation
        threshold_candidate = sorted_scores[elbow_idx]
        stable_scores = scores[scores >= threshold_candidate]
        unstable_scores = scores[scores < threshold_candidate]
        
        if len(unstable_scores) == 0 or len(stable_scores) == 0:
            return np.min(scores) - 0.01  # All stable
        
        # Check gap
        gap = np.min(stable_scores) - np.max(unstable_scores)
        if gap < 0.05:  # Very small gap
            return np.min(scores) - 0.01  # All stable
        
        threshold = threshold_candidate
    
    return threshold


def _plot_stability(stability_scores: Dict[str, float], threshold: float, 
                    component: str, species_list: List[str], 
                    stable_species: List[str], save_path: Optional[str],
                    fig_size: Tuple[float, float], font_path: Optional[str],
                    show_stats: bool, threshold_method: str):
    """
    Create bar plot visualization for stability scores.
    """
    # Set font if provided
    if font_path and Path(font_path).exists():
        font_prop = fm.FontProperties(fname=font_path)
        plt.rcParams['font.family'] = font_prop.get_name()
    
    # Create figure
    fig, ax = plt.subplots(figsize=fig_size)
    
    # Prepare data
    species_names = list(stability_scores.keys())
    scores = list(stability_scores.values())
    
    # Define colors (from compare_core_iModulon)
    stable_color = '#C1C6E8'  # Blue
    unstable_color = '#F0DDD2'  # Peach
    
    # Create bars with appropriate colors
    colors = [stable_color if species in stable_species else unstable_color 
             for species in species_names]
    
    bars = ax.bar(range(len(species_names)), scores, color=colors, edgecolor='black', linewidth=1)
    
    # Add threshold line
    ax.axhline(y=threshold, color='red', linestyle='--', linewidth=2, 
              label=f'Threshold ({threshold_method}): {threshold:.3f}')
    
    # Customize axes
    ax.set_xticks(range(len(species_names)))
    ax.set_xticklabels(species_names, rotation=45, ha='right')
    ax.set_ylabel('Stability Score (Mean Pairwise Correlation)', fontsize=12)
    ax.set_title(f'Core iModulon {component} Stability Analysis', fontsize=14)
    
    # Set y-axis limits
    ax.set_ylim([min(0, min(scores) - 0.1), 1.0])
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add legend
    legend_elements = [
        Patch(facecolor=stable_color, edgecolor='black', label=f'Stable (n={len(stable_species)})'),
        Patch(facecolor=unstable_color, edgecolor='black', label=f'Unstable (n={len(species_names)-len(stable_species)})')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Add statistics if requested
    if show_stats:
        stats_text = f"Mean: {np.mean(scores):.3f}\n"
        stats_text += f"Std: {np.std(scores):.3f}\n"
        stats_text += f"Range: [{min(scores):.3f}, {max(scores):.3f}]"
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
               verticalalignment='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Tight layout
    plt.tight_layout()
    
    # Save or show
    if save_path:
        save_path = Path(save_path)
        if save_path.is_dir():
            save_path = save_path / f"{component}_stability.svg"
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Stability plot saved to {save_path}")
    else:
        plt.show()
    
    plt.close()