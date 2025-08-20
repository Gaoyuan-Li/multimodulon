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
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)


def core_iModulon_stability(
    multimodulon,
    component: str,
    save_path: Optional[str] = None,
    fig_size: Tuple[float, float] = (6, 4),
    font_path: Optional[str] = None,
    show_stats: bool = True
) -> Tuple[List[str], float, float, Dict[str, float]]:
    """
    Quantify core iModulon stability across species using pairwise correlations.
    
    This function calculates how similar a core iModulon is across different species
    by computing the mean pairwise correlation of M matrix weights. It uses a 
    statistical approach (mean ± 2 standard deviations) to define stable ranges.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        MultiModulon instance containing the data
    component : str
        Component name (e.g., 'Core_1', 'Iron')
    save_path : str, optional
        Path to save the plot. If None, displays the plot
    fig_size : tuple, optional
        Figure size as (width, height). Default: (6, 4)
    font_path : str, optional
        Path to custom font file
    show_stats : bool, optional
        Whether to show statistics on the plot. Default: True
        
    Returns
    -------
    stable_species : list
        List of species names classified as stable (within mean ± 2*std)
    stable_min : float
        Lower boundary of stable range (mean - 2*std)
    stable_max : float  
        Upper boundary of stable range (mean + 2*std, capped at 1.0)
    stability_scores : dict
        Dictionary mapping species names to stability scores
        
    Examples
    --------
    >>> # Simple usage
    >>> stable, min_bound, max_bound, scores = mm.core_iModulon_stability("Core_1")
    
    >>> # With custom font
    >>> stable, min_bound, max_bound, scores = mm.core_iModulon_stability(
    ...     "Iron", 
    ...     font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
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
    
    # Calculate stable range using statistical approach (mean ± 2 std)
    scores_array = np.array(list(stability_scores.values()))
    mean_score = np.mean(scores_array)
    std_score = np.std(scores_array)
    
    # Define stable range boundaries
    stable_min = max(0.0, mean_score - 2 * std_score)  # Don't go below 0
    stable_max = min(1.0, mean_score + 2 * std_score)  # Don't go above 1
    
    # Classify species as stable if within the range
    stable_species = [species for species, score in stability_scores.items() 
                     if stable_min <= score <= stable_max]
    
    # Create visualization
    _plot_stability(stability_scores, stable_min, stable_max, component, species_with_component,
                   stable_species, save_path, fig_size, font_path, show_stats)
    
    return stable_species, stable_min, stable_max, stability_scores




def _plot_stability(stability_scores: Dict[str, float], stable_min: float, stable_max: float,
                    component: str, species_list: List[str], 
                    stable_species: List[str], save_path: Optional[str],
                    fig_size: Tuple[float, float], font_path: Optional[str],
                    show_stats: bool):
    """
    Create bar plot visualization for stability scores.
    """
    # Set font if provided
    font_prop = None
    if font_path and Path(font_path).exists():
        font_prop = fm.FontProperties(fname=font_path)
    
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
    
    # Add stable range lines (mean ± 2 std)
    ax.axhline(y=stable_min, color='gray', linestyle='--', linewidth=1.5, alpha=0.7,
              label=f'Stable range: {stable_min:.2f} - {stable_max:.2f}')
    ax.axhline(y=stable_max, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Fill the stable range area
    ax.axhspan(stable_min, stable_max, color='lightblue', alpha=0.1)
    
    # Customize axes
    ax.set_xticks(range(len(species_names)))
    ax.set_xticklabels(species_names, rotation=45, ha='right', fontproperties=font_prop)
    ax.set_ylabel('Stability Score (Mean Pairwise Correlation)', fontsize=12, fontproperties=font_prop)
    ax.set_title(f'Core iModulon {component} Stability Analysis', fontsize=14, fontproperties=font_prop)
    
    # Set y-axis limits
    ax.set_ylim([min(0, min(scores) - 0.1), 1.0])
    
    # Apply font to y-tick labels
    if font_prop:
        for label in ax.get_yticklabels():
            label.set_fontproperties(font_prop)
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add legend
    legend_elements = [
        Patch(facecolor=stable_color, edgecolor='black', label=f'Stable (n={len(stable_species)})'),
        Patch(facecolor=unstable_color, edgecolor='black', label=f'Unstable (n={len(species_names)-len(stable_species)})')
    ]
    legend = ax.legend(handles=legend_elements, loc='upper right')
    if font_prop:
        for text in legend.get_texts():
            text.set_fontproperties(font_prop)
    
    # Add statistics if requested
    if show_stats:
        stats_text = f"Mean: {np.mean(scores):.3f}\n"
        stats_text += f"Std: {np.std(scores):.3f}"
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
               verticalalignment='top', fontsize=10, fontproperties=font_prop,
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