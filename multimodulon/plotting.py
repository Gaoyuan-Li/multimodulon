"""Plotting functions for MultiModulon visualization."""

from __future__ import annotations

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import warnings
from typing import Dict, Optional, Tuple, List, Union
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch
import matplotlib.patches as patches
import os

# Import COG constants from core_io
from .core_io import COG_COLORS, COG_LETTER_CODES

logger = logging.getLogger(__name__)


def view_iModulon_weights(multimodulon, species: str, component: str, save_path: Optional[str] = None, 
                  fig_size: Tuple[float, float] = (6, 4), font_path: Optional[str] = None,
                  show_COG: bool = False):
    """
    Visualize gene weights for a specific iModulon component in a species.
    
    Creates a scatter plot showing gene weights across the genome, with genes
    positioned by their genomic coordinates.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        MultiModulon instance containing the data
    species : str
        Species/strain name
    component : str
        Component name (e.g., 'Core_1', 'Unique_1')
    save_path : str, optional
        Path to save the plot. Can be a directory or file path.
        If directory, saves as '{species}_{component}_iModulon.svg'
        If None, displays plot without saving
    fig_size : tuple, optional
        Figure size as (width, height). Default: (6, 4)
    font_path : str, optional
        Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
        If provided, uses this font for all text elements
    show_COG : bool, optional
        If True, color genes above threshold by their COG category. Default: False
        
    Raises
    ------
    ValueError
        If species not found, M matrix not available, or component not found
    """
    # Validate species
    if species not in multimodulon._species_data:
        raise ValueError(f"Species '{species}' not found in loaded data")
    
    species_data = multimodulon._species_data[species]
    
    # Check if M matrix exists
    if species_data.M is None:
        raise ValueError(f"M matrix not found for {species}. Please run ICA first.")
    
    # Check if component exists
    if component not in species_data.M.columns:
        raise ValueError(f"Component '{component}' not found in M matrix. "
                       f"Available components: {list(species_data.M.columns)}")
    
    # Check if gene_table exists
    if species_data.gene_table is None:
        raise ValueError(f"Gene table not found for {species}. Please run create_gene_table() first.")
    
    # Get gene weights for the component
    gene_weights = species_data.M[component]
    
    # Get gene positions from gene_table
    gene_table = species_data.gene_table
    
    # Create mapping from leftmost genes (in M matrix) to species genes (in gene_table)
    leftmost_to_species = {}
    species_to_leftmost = {}
    
    if multimodulon.combined_gene_db is not None and species in multimodulon.combined_gene_db.columns:
        # First, find the leftmost gene for each row (same logic as in alignment)
        for _, row in multimodulon.combined_gene_db.iterrows():
            # Find the leftmost (first non-None) gene in this row
            leftmost_gene = None
            for col in multimodulon.combined_gene_db.columns:
                val = row[col]
                if pd.notna(val) and val != "None" and val is not None:
                    leftmost_gene = val
                    break
            
            # Get the species-specific gene
            species_gene = row[species]
            
            if leftmost_gene and pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                leftmost_to_species[leftmost_gene] = species_gene
                species_to_leftmost[species_gene] = leftmost_gene
    
    # Create lists for plotting
    genes_for_plotting = []
    weights_for_plotting = []
    positions_for_plotting = []
    
    # Debug: log mapping statistics
    logger.debug(f"Total genes in M matrix: {len(gene_weights)}")
    logger.debug(f"Leftmost to species mappings found: {len(leftmost_to_species)}")
    
    # Process each gene in M matrix (which uses leftmost names)
    genes_not_found = []
    for gene_in_M in gene_weights.index:
        weight = gene_weights[gene_in_M]
        
        # Find the corresponding species gene name
        if gene_in_M in leftmost_to_species:
            # This is a leftmost gene, map to species gene
            species_gene = leftmost_to_species[gene_in_M]
            if species_gene in gene_table.index:
                gene_info = gene_table.loc[species_gene]
                center = (gene_info['start'] + gene_info['end']) / 2
                genes_for_plotting.append(species_gene)
                weights_for_plotting.append(weight)
                positions_for_plotting.append(center / 1e6)  # Convert to Mb
            else:
                genes_not_found.append(f"{gene_in_M}->{species_gene} (not in gene_table)")
        elif gene_in_M in gene_table.index:
            # This gene name exists directly in gene_table (not renamed)
            gene_info = gene_table.loc[gene_in_M]
            center = (gene_info['start'] + gene_info['end']) / 2
            genes_for_plotting.append(gene_in_M)
            weights_for_plotting.append(weight)
            positions_for_plotting.append(center / 1e6)  # Convert to Mb
        else:
            genes_not_found.append(f"{gene_in_M} (no mapping found)")
    
    if genes_not_found:
        logger.debug(f"Genes not plotted ({len(genes_not_found)}): {genes_not_found[:10]}...")
        logger.debug(f"Genes successfully plotted: {len(genes_for_plotting)}/{len(gene_weights)}")
    
    if not genes_for_plotting:
        raise ValueError(f"No gene position information found for genes in {species}")
    
    # Now we have the correct mapping
    x_positions = positions_for_plotting
    y_weights = weights_for_plotting
    genes_with_pos = genes_for_plotting
    
    # Get threshold if available
    threshold = None
    if hasattr(species_data, '_M_thresholds') and species_data._M_thresholds is not None:
        if component in species_data._M_thresholds.index:
            threshold = species_data._M_thresholds.loc[component, 'M_threshold']
    
    # Create figure - the given fig_size is for the plot area only
    # If showing COG, we need extra space for the legend
    if show_COG and 'COG_category' in gene_table.columns:
        # Add extra width for legend (80% more for long COG category names)
        fig_width = fig_size[0] * 1.8
        fig = plt.figure(figsize=(fig_width, fig_size[1]))
        # Create axis that uses only the original fig_size width
        ax = fig.add_axes([0.1, 0.1, fig_size[0]/fig_width * 0.8, 0.8])
    else:
        _, ax = plt.subplots(figsize=fig_size)
    
    # Set font properties if provided
    if font_path and os.path.exists(font_path):
        font_prop = fm.FontProperties(fname=font_path)
        plt.rcParams['font.family'] = font_prop.get_name()
    
    # Create scatter plot with color coding
    if show_COG and 'COG_category' in gene_table.columns and 'Description' in gene_table.columns:
        # Helper function to get COG color
        def get_cog_color(gene_name, weight):
            if threshold is not None and abs(weight) <= threshold:
                return 'grey'  # Below threshold
            
            # Try to get COG category for this gene
            cog_cat = None
            if gene_name in gene_table.index:
                cog_info = gene_table.loc[gene_name, 'COG_category']
                desc_info = gene_table.loc[gene_name, 'Description']
                
                # Check COG_category first
                if pd.notna(cog_info) and cog_info != '-':
                    # Extract first letter if it's a letter code
                    if isinstance(cog_info, str) and len(cog_info) > 0:
                        first_letter = cog_info[0]
                        if first_letter in COG_LETTER_CODES:
                            cog_cat = COG_LETTER_CODES[first_letter]
                
                # If not found in COG_category, check Description
                if cog_cat is None and pd.notna(desc_info) and desc_info != '-':
                    # Check if description matches any COG category
                    for cat_name in COG_COLORS.keys():
                        if cat_name != 'No COG annotation' and cat_name.lower() in str(desc_info).lower():
                            cog_cat = cat_name
                            break
            
            if cog_cat and cog_cat in COG_COLORS:
                return COG_COLORS[cog_cat]
            else:
                return COG_COLORS['No COG annotation']
        
        # Get colors for all genes
        colors = [get_cog_color(gene, weight) for gene, weight in zip(genes_with_pos, y_weights)]
        
        # Create scatter plot
        ax.scatter(x_positions, y_weights, alpha=0.6, s=20, c=colors)
        
        # Add horizontal threshold lines if available
        if threshold is not None:
            ax.axhline(y=threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
            ax.axhline(y=-threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
        
        # Create legend
        unique_colors = {}
        for _, color, weight in zip(genes_with_pos, colors, y_weights):
            if threshold is None or abs(weight) > threshold:
                # Get COG category name for this color
                for cat_name, cat_color in COG_COLORS.items():
                    if cat_color == color:
                        unique_colors[cat_name] = color
                        break
        
        # Sort categories by type
        info_storage = ['Translation, ribosomal structure and biogenesis', 'RNA processing and modification', 
                       'Transcription', 'Replication, recombination and repair', 'Chromatin structure and dynamics']
        cellular = ['Cell cycle control, cell division, chromosome partitioning', 'Nuclear structure', 
                   'Defense mechanisms', 'Signal transduction mechanisms', 'Cell wall/membrane/envelope biogenesis',
                   'Cell motility', 'Cytoskeleton', 'Extracellular structures', 
                   'Intracellular trafficking, secretion, and vesicular transport',
                   'Posttranslational modification, protein turnover, chaperones',
                   'Post-translational modification, protein turnover, and chaperones']
        metabolism = ['Energy production and conversion', 'Carbohydrate transport and metabolism',
                     'Amino acid transport and metabolism', 'Nucleotide transport and metabolism',
                     'Coenzyme transport and metabolism', 'Lipid transport and metabolism',
                     'Inorganic ion transport and metabolism', 
                     'Secondary metabolites biosynthesis, transport and catabolism',
                     'Secondary metabolites biosynthesis, transport, and catabolism']
        poorly_char = ['Function unknown', 'No COG annotation']
        
        sorted_categories = []
        for cat_list in [info_storage, cellular, metabolism, poorly_char]:
            for cat in cat_list:
                if cat in unique_colors:
                    sorted_categories.append(cat)
        
        if sorted_categories:
            # Create legend elements
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=unique_colors[cat], label=cat) 
                             for cat in sorted_categories]
            
            # Add legend with appropriate font - align to top edge
            if font_path and os.path.exists(font_path):
                legend = ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1),
                                 frameon=True, fontsize=10)
                for text in legend.get_texts():
                    text.set_fontproperties(font_prop)
            else:
                ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1),
                        frameon=True, fontsize=10)
    
    elif threshold is not None:
        # Original threshold-based coloring
        colors = []
        for weight in y_weights:
            if abs(weight) > threshold:
                colors.append('lightblue')  # Above threshold
            else:
                colors.append('grey')  # Below threshold
        ax.scatter(x_positions, y_weights, alpha=0.6, s=20, c=colors)
        
        # Add horizontal threshold lines
        ax.axhline(y=threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
        ax.axhline(y=-threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
    else:
        # No threshold available, use default plotting
        ax.scatter(x_positions, y_weights, alpha=0.6, s=20)
    
    # Set labels and title
    ax.set_xlabel('Gene Start (1e6)', fontsize=12)
    ax.set_ylabel('Gene Weight', fontsize=12)
    ax.set_title(f'iModulon {component} on {species}', fontsize=14)
    
    # Set font for tick labels if font_path provided
    if font_path and os.path.exists(font_path):
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontproperties(font_prop)
        ax.xaxis.label.set_fontproperties(font_prop)
        ax.yaxis.label.set_fontproperties(font_prop)
        ax.title.set_fontproperties(font_prop)
    
    # Tight layout - only apply if not showing COG (manual layout for COG)
    if not (show_COG and 'COG_category' in gene_table.columns):
        plt.tight_layout()
    
    # Save or show
    if save_path:
        # Determine save file path
        save_path = Path(save_path)
        if save_path.suffix in ['.svg', '.png', '.pdf', '.jpg']:
            # Full file path provided
            save_file = save_path
        else:
            # Directory provided, use default name
            save_path.mkdir(parents=True, exist_ok=True)
            save_file = save_path / f"{species}_{component}_iModulon.svg"
        
        # Save with tight bbox to include legend if present
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        logger.info(f"Plot saved to {save_file}")
        plt.close()
    else:
        plt.show()


def view_iModulon_activities(multimodulon, species: str, component: str, save_path: Optional[str] = None,
                            fig_size: Tuple[float, float] = (12, 3), font_path: Optional[str] = None,
                            highlight_project: Optional[Union[str, List[str]]] = None, highlight_study: Optional[str] = None,
                            highlight_condition: Optional[Union[str, List[str]]] = None,
                            show_highlight_only: bool = False,
                            show_highlight_only_color: Optional[List[str]] = None):
    """
    Visualize iModulon activities for a specific component in a species.
    
    Creates a bar plot showing component activities across samples, with samples
    grouped by project/study_accession from the sample sheet.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        MultiModulon instance containing the data
    species : str
        Species/strain name
    component : str
        Component name (e.g., 'Core_1', 'Unique_1')
    save_path : str, optional
        Path to save the plot. Can be a directory or file path.
        If directory, saves as '{species}_{component}_activities.svg'
        If None, displays plot without saving
    fig_size : tuple, optional
        Figure size as (width, height). Default: (12, 3)
    font_path : str, optional
        Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
        If provided, uses this font for all text elements
    highlight_project : str or list of str, optional
        Project name(s) to highlight with different colors. Uses 'condition' column for legend.
        Can be a single project name (string) or multiple project names (list)
    highlight_study : str, optional
        Study accession to highlight with different colors. Uses 'sample_description' column for legend
    highlight_condition : str or list of str, optional
        Condition name(s) to highlight. Only used when 'condition' column exists in sample sheet
    show_highlight_only : bool, optional
        If True, only show highlighted conditions. Default: False
    show_highlight_only_color : list of str, optional
        Colors for highlighted conditions when show_highlight_only=True.
        Must match the number of highlighted conditions
        
    Raises
    ------
    ValueError
        If species not found, A matrix not available, or component not found
    """
    # Validate species
    if species not in multimodulon._species_data:
        raise ValueError(f"Species '{species}' not found in loaded data")
    
    species_data = multimodulon._species_data[species]
    
    # Check if A matrix exists
    if species_data.A is None:
        raise ValueError(f"A matrix not found for {species}. Please run ICA first.")
    
    # Check if component exists
    if component not in species_data.A.index:
        raise ValueError(f"Component '{component}' not found in A matrix. "
                       f"Available components: {list(species_data.A.index)}")
    
    # Get component activities
    activities = species_data.A.loc[component]
    
    # Validate highlight parameters
    if show_highlight_only and not highlight_condition:
        raise ValueError("highlight_condition must be provided when show_highlight_only=True")
        
    if show_highlight_only and highlight_condition and show_highlight_only_color:
        # Convert to list if string
        if isinstance(highlight_condition, str):
            highlight_conditions = [highlight_condition]
        else:
            highlight_conditions = highlight_condition
            
        if len(show_highlight_only_color) != len(highlight_conditions):
            raise ValueError(f"Number of colors ({len(show_highlight_only_color)}) must match "
                           f"number of highlighted conditions ({len(highlight_conditions)})")
    
    # Get sample sheet
    sample_sheet = species_data.sample_sheet
    condition_mode = False
    condition_data = {}
    group_col = None
    
    if sample_sheet is None:
        # If no sample sheet, use simple bar plot without grouping
        x_positions = range(len(activities))
        x_labels = activities.index.tolist()
        group_labels = []
    else:
        # Check if condition column exists
        if 'condition' in sample_sheet.columns:
            condition_mode = True
            # Ensure activities and sample_sheet are aligned
            common_samples = [s for s in activities.index if s in sample_sheet.index]
            activities = activities[common_samples]
            sample_sheet = sample_sheet.loc[common_samples]
            
            # Group samples by condition AND project (if available)
            # This prevents grouping same conditions from different projects
            for sample in common_samples:
                condition = sample_sheet.loc[sample, 'condition']
                # Check if project column exists to create unique condition keys
                if 'project' in sample_sheet.columns:
                    project = sample_sheet.loc[sample, 'project']
                    # Create a unique key that combines condition and project
                    condition_key = (condition, project)
                else:
                    condition_key = (condition, None)
                    
                if condition_key not in condition_data:
                    condition_data[condition_key] = {
                        'samples': [],
                        'activities': [],
                        'mean_activity': 0,
                        'condition': condition,
                        'project': project if 'project' in sample_sheet.columns else None
                    }
                condition_data[condition_key]['samples'].append(sample)
                condition_data[condition_key]['activities'].append(activities[sample])
            
            # Calculate mean activities for each condition
            conditions = list(condition_data.keys())
            for condition_key in conditions:
                mean_act = np.mean(condition_data[condition_key]['activities'])
                condition_data[condition_key]['mean_activity'] = mean_act
            
            # Filter conditions if show_highlight_only is True
            if show_highlight_only and highlight_condition:
                if isinstance(highlight_condition, str):
                    highlight_conditions = [highlight_condition]
                else:
                    highlight_conditions = highlight_condition
                
                # Keep only highlighted conditions (check the condition part of the key)
                filtered_conditions = [c for c in conditions if condition_data[c]['condition'] in highlight_conditions]
                if not filtered_conditions:
                    raise ValueError(f"None of the highlighted conditions {highlight_conditions} found in data")
                conditions = filtered_conditions
            
            # Still use project/study grouping for x-axis
            group_col = None
            if 'project' in sample_sheet.columns:
                group_col = 'project'
            elif 'study_accession' in sample_sheet.columns:
                group_col = 'study_accession'
            
            if group_col:
                # Create a mapping of conditions to their project/study groups
                condition_groups = {}
                for condition_key in conditions:
                    # Use the project/study directly from the condition data
                    if group_col == 'project' and condition_data[condition_key]['project'] is not None:
                        condition_groups[condition_key] = condition_data[condition_key]['project']
                    else:
                        # Get the most common group for this condition
                        condition_samples = condition_data[condition_key]['samples']
                        groups = [sample_sheet.loc[s, group_col] for s in condition_samples if s in sample_sheet.index]
                        if groups:
                            # Use the most common group
                            most_common_group = max(set(groups), key=groups.count)
                            condition_groups[condition_key] = most_common_group
                        else:
                            condition_groups[condition_key] = 'Unknown'
                
                # Sort conditions by their project/study group and then by condition name
                sorted_conditions = sorted(conditions, key=lambda c: (condition_groups.get(c, 'Unknown'), condition_data[c]['condition']))
                conditions = sorted_conditions
                
                # Get mean activities in the sorted order
                mean_activities = [condition_data[c]['mean_activity'] for c in conditions]
                
                group_labels = [condition_groups.get(c, 'Unknown') for c in conditions]
                unique_groups = []
                for g in group_labels:
                    if g not in unique_groups:
                        unique_groups.append(g)
                
                x_positions = range(len(conditions))
                _ = [unique_groups.index(g) if g in unique_groups else -1 for g in group_labels]
            else:
                # No grouping column, just use conditions as is
                mean_activities = [condition_data[c]['mean_activity'] for c in conditions]
                x_positions = range(len(conditions))
                _ = conditions
                group_labels = []
        else:
            # Original logic for project/study grouping
            # Determine grouping column
            group_col = None
            if 'project' in sample_sheet.columns:
                group_col = 'project'
            elif 'study_accession' in sample_sheet.columns:
                group_col = 'study_accession'
            
            if group_col:
                # Ensure activities and sample_sheet are aligned
                common_samples = [s for s in activities.index if s in sample_sheet.index]
                activities = activities[common_samples]
                sample_sheet = sample_sheet.loc[common_samples]
                
                # Get group labels
                group_labels = sample_sheet[group_col].fillna('Unknown').astype(str).tolist()
                unique_groups = []
                for g in group_labels:
                    if g not in unique_groups:
                        unique_groups.append(g)
                
                x_positions = range(len(activities))
                _ = [unique_groups.index(g) if g in unique_groups else -1 for g in group_labels]
            else:
                # No grouping column available
                x_positions = range(len(activities))
                _ = activities.index.tolist()
                group_labels = []
    
    # Prepare colors for bars
    colors = []
    legend_elements = []
    
    if condition_mode and show_highlight_only and show_highlight_only_color:
        # Use provided colors for highlighted conditions
        for i, condition_key in enumerate(conditions):
            color = show_highlight_only_color[i]
            colors.append(color)
            # Get the condition name for legend
            condition_name = condition_data[condition_key]['condition']
            legend_elements.append((condition_name, color))
    elif condition_mode and highlight_condition and not show_highlight_only:
        # Highlight specific conditions but show all
        if isinstance(highlight_condition, str):
            highlight_conditions = [highlight_condition]
        else:
            highlight_conditions = highlight_condition
            
        # Use tab20 for more colors (20 distinct colors)
        colormap = plt.cm.tab20
        color_idx = 0
        for condition_key in conditions:
            condition_name = condition_data[condition_key]['condition']
            if condition_name in highlight_conditions:
                # Check if we need to add project suffix
                project = condition_data[condition_key]['project']
                if project:
                    legend_label = f"{condition_name} ({project})"
                else:
                    legend_label = condition_name
                    
                if legend_label not in [elem[0] for elem in legend_elements]:
                    color = colormap(color_idx % 20)
                    legend_elements.append((legend_label, color))
                    colors.append(color)
                    color_idx += 1
                else:
                    # Find existing color for this condition
                    for label, col in legend_elements:
                        if label == legend_label:
                            colors.append(col)
                            break
            else:
                colors.append('lightblue')
    elif condition_mode and highlight_project and sample_sheet is not None and 'project' in sample_sheet.columns:
        # In condition mode with highlighting, color by condition
        # Convert single string to list for uniform handling
        if isinstance(highlight_project, str):
            highlight_projects = [highlight_project]
        else:
            highlight_projects = highlight_project
        
        for i, condition_key in enumerate(conditions):
            condition_name = condition_data[condition_key]['condition']
            project = condition_data[condition_key]['project']
            
            # Check if this condition's project is highlighted
            condition_highlighted = project and any(str(hp) == str(project) for hp in highlight_projects)
            
            if condition_highlighted:
                # Create legend label with project name if there might be duplicates
                # Check if this condition appears in multiple highlighted projects
                condition_appears_multiple = False
                for other_key in conditions:
                    if (condition_data[other_key]['condition'] == condition_name and 
                        other_key != condition_key and
                        condition_data[other_key]['project'] and
                        any(str(hp) == str(condition_data[other_key]['project']) for hp in highlight_projects)):
                        condition_appears_multiple = True
                        break
                
                if condition_appears_multiple:
                    legend_label = f"{condition_name} ({project})"
                else:
                    legend_label = condition_name
                
                # Assign unique color per condition using tab20 for more colors
                if legend_label not in [elem[0] for elem in legend_elements]:
                    color = plt.cm.tab20(len(legend_elements) % 20)
                    legend_elements.append((legend_label, color))
                    colors.append(color)
                else:
                    # Find existing color for this condition
                    for label, col in legend_elements:
                        if label == legend_label:
                            colors.append(col)
                            break
            else:
                colors.append('lightblue')
    elif condition_mode:
        # In condition mode without highlighting, use default colors
        colors = ['lightblue'] * len(conditions)
    elif highlight_project and sample_sheet is not None and 'project' in sample_sheet.columns:
        # Convert single string to list for uniform handling
        if isinstance(highlight_project, str):
            highlight_projects = [highlight_project]
        else:
            highlight_projects = highlight_project
        
        # Color by project with condition as legend
        for i, sample in enumerate(activities.index):
            if sample in sample_sheet.index:
                project = str(sample_sheet.loc[sample, 'project'])
                # Check if project is in the highlight list
                project_highlighted = any(str(hp) == project for hp in highlight_projects)
                
                if project_highlighted:
                    condition = sample_sheet.loc[sample, 'condition'] if 'condition' in sample_sheet.columns else 'Unknown'
                    # Assign unique color per condition using tab20 for more colors
                    if condition not in [elem[0] for elem in legend_elements]:
                        color = plt.cm.tab20(len(legend_elements) % 20)
                        legend_elements.append((condition, color))
                        colors.append(color)
                    else:
                        # Find existing color for this condition
                        for cond, col in legend_elements:
                            if cond == condition:
                                colors.append(col)
                                break
                else:
                    colors.append('lightblue')
            else:
                colors.append('lightblue')
    elif highlight_study and sample_sheet is not None and 'study_accession' in sample_sheet.columns:
        # Color by study with sample_description as legend
        for i, sample in enumerate(activities.index):
            if sample in sample_sheet.index:
                study = sample_sheet.loc[sample, 'study_accession']
                if str(study) == str(highlight_study):
                    desc = sample_sheet.loc[sample, 'sample_description'] if 'sample_description' in sample_sheet.columns else 'Unknown'
                    # Assign unique color per description using tab20 for more colors
                    if desc not in [elem[0] for elem in legend_elements]:
                        color = plt.cm.tab20(len(legend_elements) % 20)
                        legend_elements.append((desc, color))
                        colors.append(color)
                    else:
                        # Find existing color for this description
                        for d, col in legend_elements:
                            if d == desc:
                                colors.append(col)
                                break
                else:
                    colors.append('lightblue')
            else:
                colors.append('lightblue')
    else:
        # Default light blue color for all bars
        colors = ['lightblue'] * len(activities)
    
    # Create figure
    _, ax = plt.subplots(figsize=fig_size)
    
    # Set font properties if provided
    if font_path and os.path.exists(font_path):
        font_prop = fm.FontProperties(fname=font_path)
        plt.rcParams['font.family'] = font_prop.get_name()
    
    # Create bar plot
    if condition_mode:
        # Plot averaged bars for conditions with reduced gap between bars
        ax.bar(x_positions, mean_activities, width=0.75, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # Add individual sample points as black dots
        for i, condition_key in enumerate(conditions):
            sample_activities = condition_data[condition_key]['activities']
            # All dots at the horizontal center of the bar
            x_points = [i] * len(sample_activities)
            ax.scatter(x_points, sample_activities, color='black', s=10, zorder=10, alpha=0.7)
    else:
        # Original bar plot for individual samples
        ax.bar(x_positions, activities.values, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Add horizontal line at y=0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    
    # Reduce gap between bars by adjusting x-axis limits
    if condition_mode:
        # Set tighter x-axis limits to reduce gaps
        ax.set_xlim(-0.5, len(conditions) - 0.5)
    
    # Remove grid
    ax.grid(False)
    
    # Add vertical dotted lines to separate projects
    if group_col and group_labels:
        current_group = group_labels[0]
        for i in range(1, len(group_labels)):
            if group_labels[i] != current_group:
                # Add vertical line at boundary
                ax.axvline(x=i-0.5, color='gray', linestyle=':', linewidth=1, alpha=0.7)
                current_group = group_labels[i]
    
    # Set labels and title
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('iModulon Activity', fontsize=12)
    ax.set_title(f'iModulon {component} on {species}', fontsize=14)
    
    # Set x-axis ticks
    if condition_mode and show_highlight_only:
        # No x-axis ticks when showing highlight only
        ax.set_xticks([])
    elif condition_mode and group_col and group_labels:
        # In condition mode with grouping, show group labels like normal
        # Show group labels on x-axis
        tick_positions = []
        tick_labels = []
        current_group = None
        group_start = 0
        
        for i, group in enumerate(group_labels):
            if group != current_group:
                if current_group is not None:
                    # Add tick for previous group
                    tick_positions.append((group_start + i - 1) / 2)
                    tick_labels.append(current_group)
                current_group = group
                group_start = i
        
        # Add last group
        if current_group is not None:
            tick_positions.append((group_start + len(group_labels) - 1) / 2)
            tick_labels.append(current_group)
        
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=45, ha='right')
    elif condition_mode:
        # In condition mode without grouping, don't show ticks
        ax.set_xticks([])
    elif group_col and group_labels:
        # Show group labels on x-axis
        tick_positions = []
        tick_labels = []
        current_group = None
        group_start = 0
        
        for i, group in enumerate(group_labels):
            if group != current_group:
                if current_group is not None:
                    # Add tick for previous group
                    tick_positions.append((group_start + i - 1) / 2)
                    tick_labels.append(current_group)
                current_group = group
                group_start = i
        
        # Add last group
        if current_group is not None:
            tick_positions.append((group_start + len(group_labels) - 1) / 2)
            tick_labels.append(current_group)
        
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=45, ha='right')
    else:
        # No grouping - don't show individual sample names
        ax.set_xticks([])
    
    # Add legend if highlighting
    if legend_elements:
        from matplotlib.patches import Patch
        patches = [Patch(color=color, label=label) for label, color in legend_elements]
        
        # Calculate number of columns needed based on figure height
        # Base: 8 items per column for fig_height=3
        n_items = len(legend_elements)
        fig_height = fig_size[1]  # Get the figure height
        items_per_column = round(8 * (fig_height / 3))  # Scale proportionally
        
        # Calculate number of columns needed
        ncol = max(1, (n_items + items_per_column - 1) // items_per_column)  # Ceiling division
        
        # Create legend with multiple columns
        legend = ax.legend(handles=patches, loc='center left', bbox_to_anchor=(1.02, 0.5),
                          ncol=ncol, frameon=True, fontsize=9)
        
        # Apply font to legend if provided
        if font_path and os.path.exists(font_path):
            for text in legend.get_texts():
                text.set_fontproperties(font_prop)
    
    # Set font for tick labels if font_path provided
    if font_path and os.path.exists(font_path):
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontproperties(font_prop)
        ax.xaxis.label.set_fontproperties(font_prop)
        ax.yaxis.label.set_fontproperties(font_prop)
        ax.title.set_fontproperties(font_prop)
    
    # Adjust layout based on legend
    if legend_elements:
        # Calculate space needed for legend based on number of columns
        # More columns need more horizontal space
        # Base formula: right_space = 0.85 - (ncol - 1) * 0.1
        # But don't go below 0.45
        right_space = max(0.45, 0.85 - (ncol - 1) * 0.1)
        
        # Use subplots_adjust instead of tight_layout to avoid warnings
        plt.subplots_adjust(left=0.1, right=right_space, top=0.95, bottom=0.15)
    else:
        plt.tight_layout()
    
    # Save or show
    if save_path:
        # Determine save file path
        save_path = Path(save_path)
        if save_path.suffix in ['.svg', '.png', '.pdf', '.jpg']:
            # Full file path provided
            save_file = save_path
        else:
            # Directory provided, use default name
            save_path.mkdir(parents=True, exist_ok=True)
            save_file = save_path / f"{species}_{component}_activities.svg"
        
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        logger.info(f"Plot saved to {save_file}")
        plt.close()
    else:
        plt.show()


def view_iModulon_genes(multimodulon, species: str, component: str) -> pd.DataFrame:
    """
    Return a subset of the gene table for genes in a specific iModulon component.
    
    Uses the presence_matrix to identify which genes belong to the specified component,
    then returns the corresponding rows from the gene_table.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        MultiModulon instance containing the data
    species : str
        Species/strain name
    component : str
        Component name (e.g., 'Core_1', 'Unique_1')
        
    Returns
    -------
    pd.DataFrame
        Subset of gene_table containing only genes in the specified component.
        Returns empty DataFrame if no genes found or if data not available.
        
    Raises
    ------
    ValueError
        If species not found or if presence_matrix/gene_table not available
        
    Examples
    --------
    >>> gene_subset = view_iModulon_genes(mm, 'E_coli', 'Core_1')
    >>> print(f"Found {len(gene_subset)} genes in Core_1")
    """
    # Validate species
    if species not in multimodulon._species_data:
        raise ValueError(f"Species '{species}' not found in loaded data")
    
    species_data = multimodulon._species_data[species]
    
    # Check if presence_matrix exists
    if species_data._presence_matrix is None:
        raise ValueError(f"presence_matrix not found for {species}. Please run optimize_M_thresholds() first.")
    
    # Check if gene_table exists
    if species_data.gene_table is None:
        raise ValueError(f"gene_table not found for {species}.")
    
    # Check if component exists
    if component not in species_data._presence_matrix.columns:
        raise ValueError(f"Component '{component}' not found in presence_matrix. "
                       f"Available components: {list(species_data._presence_matrix.columns)}")
    
    # Get genes in this component from presence_matrix
    component_genes = species_data._presence_matrix[
        species_data._presence_matrix[component] == 1
    ].index.tolist()
    
    if not component_genes:
        logger.warning(f"No genes found in component '{component}' for species '{species}'")
        return pd.DataFrame()
    
    # Map component genes to gene_table
    # First check if we need to map through combined_gene_db
    genes_in_table = []
    
    if multimodulon.combined_gene_db is not None and species in multimodulon.combined_gene_db.columns:
        # Map from leftmost genes to species-specific genes
        for gene in component_genes:
            # Check if gene exists directly in gene_table first
            if gene in species_data.gene_table.index:
                genes_in_table.append(gene)
            else:
                # Try to map through combined_gene_db
                for _, row in multimodulon.combined_gene_db.iterrows():
                    # Find leftmost gene
                    leftmost_gene = None
                    for col in multimodulon.combined_gene_db.columns:
                        val = row[col]
                        if pd.notna(val) and val != "None" and val is not None:
                            leftmost_gene = val
                            break
                    
                    # If this is our gene, get the species-specific name
                    if leftmost_gene == gene:
                        species_gene = row[species]
                        if pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                            if species_gene in species_data.gene_table.index:
                                genes_in_table.append(species_gene)
                                break
    else:
        # No combined_gene_db, just check directly
        genes_in_table = [g for g in component_genes if g in species_data.gene_table.index]
    
    # Get subset of gene_table
    if genes_in_table:
        gene_subset = species_data.gene_table.loc[genes_in_table].copy()
        
        # Sort by start position if the column exists
        if 'start' in gene_subset.columns:
            gene_subset = gene_subset.sort_values('start')
            
        logger.info(f"Found {len(gene_subset)} genes in component '{component}' for species '{species}'")
        return gene_subset
    else:
        logger.warning(f"No genes from component '{component}' found in gene_table for species '{species}'")
        return pd.DataFrame()


def view_core_iModulon_weights(multimodulon, component: str, save_path: Optional[str] = None,
                       fig_size: Tuple[float, float] = (6, 4), font_path: Optional[str] = None,
                       show_COG: bool = False, reference_order: Optional[List[str]] = None):
    """
    Visualize a core iModulon component across all species.
    
    Creates individual plots for each species showing the same core component,
    and saves them to the specified directory. This method calls view_iModulon_weights
    for each species with the given core component.
    
    Parameters
    ----------
    multimodulon : MultiModulon
        MultiModulon instance containing the data
    component : str
        Core component name (e.g., 'Core_1', 'Core_2')
    save_path : str, optional
        Directory path to save all plots. If None, displays plots without saving.
        Each plot is saved as '{species}_{component}_iModulon.svg'
    fig_size : tuple, optional
        Figure size for individual plots as (width, height). Default: (6, 4)
    font_path : str, optional
        Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
        If provided, uses this font for all text elements in all plots
    show_COG : bool, optional
        If True, color genes above threshold by their COG category. Default: False
    reference_order : list of str, optional
        Custom order for arranging species in subplots when show_COG is True.
        First 3 species will be placed in the first row, remaining in the second row.
        Example: ['MG1655', 'BL21', 'C', 'Crooks', 'W', 'W3110']
        If not provided, species are plotted in their default order.
        
    Raises
    ------
    ValueError
        If component is not a core component or not found in any species
    """
    # Get all species
    species_list = list(multimodulon._species_data.keys())
    
    # Verify this is a core component by checking it exists in multiple species
    species_with_component = []
    for species in species_list:
        species_data = multimodulon._species_data[species]
        if species_data.M is not None and component in species_data.M.columns:
            species_with_component.append(species)
    
    if not species_with_component:
        raise ValueError(f"Component '{component}' not found in any species")
    
    if not component.startswith('Core_'):
        warnings.warn(f"Component '{component}' does not appear to be a core component "
                     f"(core components typically start with 'Core_')")
    
    # Create save directory if needed
    if save_path:
        save_path = Path(save_path)
        if not save_path.suffix:  # It's a directory
            save_path.mkdir(parents=True, exist_ok=True)
    
    # Generate plots for each species
    logger.info(f"Generating plots for core component '{component}' across {len(species_with_component)} species")
    
    if show_COG:
        # When showing COG, create a combined plot with subplots and single legend
        
        # Apply reference_order if provided
        if reference_order:
            # Filter reference_order to only include species that have the component
            ordered_species = [sp for sp in reference_order if sp in species_with_component]
            # Add any remaining species not in reference_order
            remaining_species = [sp for sp in species_with_component if sp not in ordered_species]
            species_with_component = ordered_species + remaining_species
        
        n_species = len(species_with_component)
        
        # Determine subplot layout
        if reference_order and n_species > 3:
            # Custom layout: first 3 in first row, rest in second row
            n_rows = 2
            n_cols = max(3, n_species - 3)
        elif n_species <= 3:
            n_rows = 1
            n_cols = n_species
        elif n_species <= 6:
            n_rows = 2
            n_cols = 3
        else:
            n_cols = 3
            n_rows = (n_species + 2) // 3
        
        # Create figure with subplots
        fig_width = fig_size[0] * n_cols
        fig_height = fig_size[1] * n_rows + 2  # Extra space for legend
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
        
        # Ensure axes is always a 2D array
        if n_species == 1:
            axes = np.array([[axes]])
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_cols == 1:
            axes = axes.reshape(-1, 1)
        
        # Set font properties if provided
        if font_path and os.path.exists(font_path):
            font_prop = fm.FontProperties(fname=font_path)
            plt.rcParams['font.family'] = font_prop.get_name()
        
        # Track all unique COG categories across all species
        all_unique_colors = {}
        
        # Plot each species
        for idx, species in enumerate(species_with_component):
            if reference_order and n_species > 3:
                # Custom layout for reference_order
                if idx < 3:
                    # First 3 species in first row
                    row = 0
                    col = idx
                else:
                    # Remaining species in second row
                    row = 1
                    col = idx - 3
            else:
                # Default layout
                row = idx // n_cols
                col = idx % n_cols
            ax = axes[row, col]
            
            species_data = multimodulon._species_data[species]
            
            # Get gene weights and gene table
            gene_weights = species_data.M[component]
            gene_table = species_data.gene_table
            
            # Create mapping from leftmost genes (in M matrix) to species genes (in gene_table)
            leftmost_to_species = {}
            
            if multimodulon.combined_gene_db is not None and species in multimodulon.combined_gene_db.columns:
                # First, find the leftmost gene for each row (same logic as in alignment)
                for _, row in multimodulon.combined_gene_db.iterrows():
                    # Find the leftmost (first non-None) gene in this row
                    leftmost_gene = None
                    for col in multimodulon.combined_gene_db.columns:
                        val = row[col]
                        if pd.notna(val) and val != "None" and val is not None:
                            leftmost_gene = val
                            break
                    
                    # Get the species-specific gene
                    species_gene = row[species]
                    
                    if leftmost_gene and pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                        leftmost_to_species[leftmost_gene] = species_gene
            
            # Create lists for plotting
            genes_for_plotting = []
            weights_for_plotting = []
            positions_for_plotting = []
            
            # Process each gene in M matrix (which uses leftmost names)
            for gene_in_M in gene_weights.index:
                weight = gene_weights[gene_in_M]
                
                # Find the corresponding species gene name
                if gene_in_M in leftmost_to_species:
                    # This is a leftmost gene, map to species gene
                    species_gene = leftmost_to_species[gene_in_M]
                    if species_gene in gene_table.index:
                        gene_info = gene_table.loc[species_gene]
                        center = (gene_info['start'] + gene_info['end']) / 2
                        genes_for_plotting.append(species_gene)
                        weights_for_plotting.append(weight)
                        positions_for_plotting.append(center / 1e6)  # Convert to Mb
                elif gene_in_M in gene_table.index:
                    # This gene name exists directly in gene_table (not renamed)
                    gene_info = gene_table.loc[gene_in_M]
                    center = (gene_info['start'] + gene_info['end']) / 2
                    genes_for_plotting.append(gene_in_M)
                    weights_for_plotting.append(weight)
                    positions_for_plotting.append(center / 1e6)  # Convert to Mb
            
            # Now we have the correct mapping
            x_positions = positions_for_plotting
            y_weights = weights_for_plotting
            genes_with_pos = genes_for_plotting
            
            # Get threshold
            threshold = None
            if hasattr(species_data, '_M_thresholds') and species_data._M_thresholds is not None:
                if component in species_data._M_thresholds.index:
                    threshold = species_data._M_thresholds.loc[component, 'M_threshold']
            
            # Helper function to get COG color
            def get_cog_color(gene_name, weight):
                if threshold is not None and abs(weight) <= threshold:
                    return 'grey'
                
                cog_cat = None
                if gene_name in gene_table.index and 'COG_category' in gene_table.columns:
                    cog_info = gene_table.loc[gene_name, 'COG_category']
                    desc_info = gene_table.loc[gene_name, 'Description'] if 'Description' in gene_table.columns else None
                    
                    if pd.notna(cog_info) and cog_info != '-':
                        if isinstance(cog_info, str) and len(cog_info) > 0:
                            first_letter = cog_info[0]
                            if first_letter in COG_LETTER_CODES:
                                cog_cat = COG_LETTER_CODES[first_letter]
                    
                    if cog_cat is None and desc_info and pd.notna(desc_info) and desc_info != '-':
                        for cat_name in COG_COLORS.keys():
                            if cat_name != 'No COG annotation' and cat_name.lower() in str(desc_info).lower():
                                cog_cat = cat_name
                                break
                
                if cog_cat and cog_cat in COG_COLORS:
                    return COG_COLORS[cog_cat]
                else:
                    return COG_COLORS['No COG annotation']
            
            # Get colors
            colors = [get_cog_color(gene, weight) for gene, weight in zip(genes_with_pos, y_weights)]
            
            # Track unique colors for legend
            for _, color, weight in zip(genes_with_pos, colors, y_weights):
                if threshold is None or abs(weight) > threshold:
                    for cat_name, cat_color in COG_COLORS.items():
                        if cat_color == color:
                            all_unique_colors[cat_name] = color
                            break
            
            # Create scatter plot
            ax.scatter(x_positions, y_weights, alpha=0.6, s=20, c=colors)
            
            # Add threshold lines
            if threshold is not None:
                ax.axhline(y=threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
                ax.axhline(y=-threshold, color='black', linestyle=':', linewidth=1, alpha=0.7)
            
            # Set labels and title
            ax.set_xlabel('Gene Start (1e6)', fontsize=10)
            ax.set_ylabel('Gene Weight', fontsize=10)
            ax.set_title(f'{species}', fontsize=12)
            
            # Set font for labels if provided
            if font_path and os.path.exists(font_path):
                for label in ax.get_xticklabels() + ax.get_yticklabels():
                    label.set_fontproperties(font_prop)
                ax.xaxis.label.set_fontproperties(font_prop)
                ax.yaxis.label.set_fontproperties(font_prop)
                ax.title.set_fontproperties(font_prop)
        
        # Hide empty subplots
        for idx in range(n_species, n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].axis('off')
        
        # Add main title with reduced gap
        fig.suptitle(f'Core iModulon {component}', fontsize=16, y=0.99)
        if font_path and os.path.exists(font_path):
            fig.suptitle(f'Core iModulon {component}', fontsize=16, y=0.99, fontproperties=font_prop)
        
        # Create legend at bottom
        if all_unique_colors:
            # Sort categories
            info_storage = ['Translation, ribosomal structure and biogenesis', 'RNA processing and modification', 
                           'Transcription', 'Replication, recombination and repair', 'Chromatin structure and dynamics']
            cellular = ['Cell cycle control, cell division, chromosome partitioning', 'Nuclear structure', 
                       'Defense mechanisms', 'Signal transduction mechanisms', 'Cell wall/membrane/envelope biogenesis',
                       'Cell motility', 'Cytoskeleton', 'Extracellular structures', 
                       'Intracellular trafficking, secretion, and vesicular transport',
                       'Posttranslational modification, protein turnover, chaperones',
                       'Post-translational modification, protein turnover, and chaperones']
            metabolism = ['Energy production and conversion', 'Carbohydrate transport and metabolism',
                         'Amino acid transport and metabolism', 'Nucleotide transport and metabolism',
                         'Coenzyme transport and metabolism', 'Lipid transport and metabolism',
                         'Inorganic ion transport and metabolism', 
                         'Secondary metabolites biosynthesis, transport and catabolism',
                         'Secondary metabolites biosynthesis, transport, and catabolism']
            poorly_char = ['Function unknown', 'No COG annotation']
            
            sorted_categories = []
            for cat_list in [info_storage, cellular, metabolism, poorly_char]:
                for cat in cat_list:
                    if cat in all_unique_colors:
                        sorted_categories.append(cat)
            
            # Create legend elements
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=all_unique_colors[cat], label=cat) 
                             for cat in sorted_categories]
            
            # Add legend at bottom
            if font_path and os.path.exists(font_path):
                legend = fig.legend(handles=legend_elements, loc='lower center', 
                                  bbox_to_anchor=(0.5, -0.03), ncol=3, frameon=True, fontsize=10)
                for text in legend.get_texts():
                    text.set_fontproperties(font_prop)
            else:
                fig.legend(handles=legend_elements, loc='lower center', 
                         bbox_to_anchor=(0.5, -0.03), ncol=3, frameon=True, fontsize=10)
        
        # Adjust layout with reduced top margin
        plt.tight_layout(rect=[0, 0.08, 1, 0.97])
        
        # Save or show
        if save_path:
            save_path = Path(save_path)
            if save_path.suffix in ['.svg', '.png', '.pdf', '.jpg']:
                save_file = save_path
            else:
                save_path.mkdir(parents=True, exist_ok=True)
                save_file = save_path / f"Core_{component}_all_species_COG.svg"
            
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            logger.info(f"Combined plot saved to {save_file}")
            plt.close()
        else:
            plt.show()
    
    else:
        # Original behavior: individual plots
        for species in species_with_component:
            try:
                view_iModulon_weights(
                    multimodulon=multimodulon,
                    species=species,
                    component=component,
                    save_path=save_path,
                    fig_size=fig_size,
                    font_path=font_path,
                    show_COG=show_COG
                )
                logger.info(f"✓ Generated plot for {species}")
            except Exception as e:
                logger.warning(f"Failed to generate plot for {species}: {str(e)}")
    
    logger.info(f"Completed generating plots for core component '{component}'")


def compare_core_iModulon(multimodulon, component: str, y_label: str = 'Species', 
                          save_path: Optional[str] = None, fig_size: Tuple[float, float] = (12, 8),
                          font_path: Optional[str] = None, reference_order: Optional[List[str]] = None) -> Dict[str, List[str]]:
    """
    Compare a core iModulon component across species with a dual-layer heatmap.
    
    Creates a heatmap showing:
    - Component membership (filled blocks from presence_matrix)
    - Gene presence across species (edge visibility from combined_gene_db)
    
    The heatmap is divided into two sections:
    - Left: genes shared by all species
    - Right: genes not shared by all species (with green edges showing presence)
    
    Parameters
    ----------
    multimodulon : MultiModulon
        MultiModulon instance containing the data
    component : str
        Core component name (e.g., 'Core_1', 'Core_2')
    y_label : str, optional
        Y-axis label. Default: 'Species'
    save_path : str, optional
        Path to save the plot. If None, displays the plot
    fig_size : tuple, optional
        Figure size as (width, height). Default: (12, 8)
    font_path : str, optional
        Path to font file (e.g., '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf')
        If provided, uses this font for all text elements
    reference_order : list of str, optional
        Custom order for species on y-axis (top to bottom).
        Example: ['MG1655', 'W', 'W3110', 'BL21', 'C', 'Crooks']
        If not provided, species are displayed in their default order.
        
    Returns
    -------
    Dict[str, List[str]]
        Dictionary with species names as keys and lists of gene names (M matrix indices) as values
        
    Raises
    ------
    ValueError
        If component is not found in any species
    """
    # Get all species with this component
    species_with_component = []
    result_dict = {}
    
    for species in multimodulon._species_data.keys():
        species_data = multimodulon._species_data[species]
        if species_data.M is not None and component in species_data.M.columns:
            species_with_component.append(species)
            
            # Get genes in this component for this species (using presence_matrix)
            if species_data._presence_matrix is not None:
                component_genes = species_data._presence_matrix[
                    species_data._presence_matrix[component] == 1
                ].index.tolist()
                result_dict[species] = component_genes
    
    if not species_with_component:
        raise ValueError(f"Component '{component}' not found in any species")
    
    if not component.startswith('Core_'):
        warnings.warn(f"Component '{component}' does not appear to be a core component")
    
    # Apply reference_order if provided
    if reference_order:
        # Filter to only include species that have the component
        ordered_species = [sp for sp in reference_order if sp in species_with_component]
        # Add any remaining species not in reference_order
        remaining_species = [sp for sp in species_with_component if sp not in ordered_species]
        species_with_component = ordered_species + remaining_species
    
    # Set font properties if provided
    if font_path and os.path.exists(font_path):
        font_prop = fm.FontProperties(fname=font_path)
        plt.rcParams['font.family'] = font_prop.get_name()
    
    # Collect all unique genes across species for this component
    all_genes = set()
    for species in species_with_component:
        if species in result_dict:
            all_genes.update(result_dict[species])
    all_genes = sorted(list(all_genes))
    
    # Build gene presence information across species
    gene_species_presence = {}  # gene -> set of species that have this gene
    for gene in all_genes:
        gene_species_presence[gene] = set()
        
        if multimodulon.combined_gene_db is not None:
            # Check each row in combined_gene_db
            for _, row in multimodulon.combined_gene_db.iterrows():
                # Find leftmost gene
                leftmost_gene = None
                for col in multimodulon.combined_gene_db.columns:
                    val = row[col]
                    if pd.notna(val) and val != "None" and val is not None:
                        leftmost_gene = val
                        break
                
                # If this row contains our gene
                if leftmost_gene == gene:
                    # Check which species have this gene
                    for species in species_with_component:
                        if species in multimodulon.combined_gene_db.columns:
                            species_gene = row[species]
                            if pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                                gene_species_presence[gene].add(species)
                    break
    
    # Separate genes into two groups
    genes_in_all = []
    genes_not_in_all = []
    
    for gene in all_genes:
        if len(gene_species_presence.get(gene, set())) == len(species_with_component):
            genes_in_all.append(gene)
        else:
            genes_not_in_all.append(gene)
    
    # Combine with gap
    gap_width = 1
    ordered_genes = genes_in_all + genes_not_in_all
    gap_position = len(genes_in_all) if genes_in_all else -1
    
    # Create figure
    fig, ax = plt.subplots(figsize=fig_size)
    
    # Prepare data for visualization
    n_species = len(species_with_component)
    n_genes = len(ordered_genes)
    total_width = n_genes + (gap_width if gap_position > 0 else 0)
    
    # Function to create merged rectangles for adjacent species
    def create_merged_rectangles(_, x_pos, species_list, in_component_flags, exists_flags, is_right_section=False):
        rectangles = []
        i = 0
        while i < len(species_list):
            # For right section, we need to handle non-existent genes too
            if not exists_flags[i] and is_right_section:
                # Gene doesn't exist in this species - add white rectangle for right section
                start_i = i
                # Find consecutive species where gene doesn't exist
                while i < len(species_list) and not exists_flags[i]:
                    i += 1
                
                height = i - start_i
                y_pos = len(species_list) - i  # Flip to have first species at top
                
                # White rectangle for non-existent genes
                rect = patches.Rectangle((x_pos, y_pos), 1, height, linewidth=0,
                                       edgecolor='none', facecolor='white')
                rectangles.append(rect)
            elif not exists_flags[i]:
                # Left section - skip non-existent genes
                i += 1
                continue
            else:
                # Gene exists - process as before
                # Start of a group
                start_i = i
                start_in_component = in_component_flags[i]
                
                # Find consecutive species with same status
                while i < len(species_list) and exists_flags[i] and in_component_flags[i] == start_in_component:
                    i += 1
                
                # Create merged rectangle
                height = i - start_i
                y_pos = len(species_list) - i  # Flip to have first species at top
                
                if start_in_component:
                    # In component - filled light blue, no edge
                    rect = patches.Rectangle((x_pos, y_pos), 1, height, linewidth=0,
                                           edgecolor='none', facecolor='#CCE5FF')  # Light blue (slightly deeper)
                else:
                    # Not in component but exists - very light grey fill, no edge
                    rect = patches.Rectangle((x_pos, y_pos), 1, height, linewidth=0,
                                           edgecolor='none', facecolor='#F0F0F0')  # Very light grey
                rectangles.append(rect)
        
        return rectangles
    
    # Create patches for the heatmap
    all_rectangles = []
    
    # Process each gene
    for j, gene in enumerate(ordered_genes):
        # Calculate x position accounting for gap
        if gap_position > 0 and j >= gap_position:
            x_pos = j + gap_width
        else:
            x_pos = j
        
        # Determine if this is in the right section
        is_right_section = gene in genes_not_in_all
        
        # Collect status for all species
        in_component_flags = []
        exists_flags = []
        
        for species in species_with_component:
            in_component = gene in result_dict.get(species, [])
            exists_in_species = species in gene_species_presence.get(gene, set())
            
            in_component_flags.append(in_component)
            exists_flags.append(exists_in_species)
        
        # Create merged rectangles
        rectangles = create_merged_rectangles(j, x_pos, species_with_component, 
                                             in_component_flags, exists_flags, is_right_section)
        all_rectangles.extend(rectangles)
    
    # Add all rectangles to the plot
    for rect in all_rectangles:
        ax.add_patch(rect)
    
    # Set axis properties with small margin to show box edges
    margin = 0.1
    ax.set_xlim(-margin, total_width + margin)
    ax.set_ylim(-margin, n_species + margin)
    ax.set_aspect('equal')
    
    # No vertical line needed - gap will separate the sections
    
    # Set ticks and labels
    x_positions = []
    for j in range(len(ordered_genes)):
        if gap_position > 0 and j >= gap_position:
            x_positions.append(j + gap_width + 0.5)
        else:
            x_positions.append(j + 0.5)
    
    ax.set_xticks(x_positions)
    ax.set_yticks(np.arange(n_species) + 0.5)
    
    # Get gene names for x-axis labels
    x_labels = []
    for gene in ordered_genes:
        # Try to get gene_name from gene_table
        gene_name = gene  # Default to using the index
        for species in species_with_component:
            species_data = multimodulon._species_data[species]
            if species_data.gene_table is not None:
                # Map from M matrix index to species gene
                if multimodulon.combined_gene_db is not None and species in multimodulon.combined_gene_db.columns:
                    # Find the species-specific gene name
                    for _, row in multimodulon.combined_gene_db.iterrows():
                        leftmost_gene = None
                        for col in multimodulon.combined_gene_db.columns:
                            val = row[col]
                            if pd.notna(val) and val != "None" and val is not None:
                                leftmost_gene = val
                                break
                        
                        if leftmost_gene == gene:
                            species_gene = row[species]
                            if pd.notna(species_gene) and species_gene != "None" and species_gene is not None:
                                if species_gene in species_data.gene_table.index:
                                    if 'gene_name' in species_data.gene_table.columns:
                                        name = species_data.gene_table.loc[species_gene, 'gene_name']
                                        if pd.notna(name) and name != '':
                                            gene_name = name
                                            break
                
                # Also check if gene exists directly in gene_table
                elif gene in species_data.gene_table.index:
                    if 'gene_name' in species_data.gene_table.columns:
                        name = species_data.gene_table.loc[gene, 'gene_name']
                        if pd.notna(name) and name != '':
                            gene_name = name
                            break
        
        x_labels.append(gene_name)
    
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    # Reverse y-axis labels to match top-to-bottom order
    ax.set_yticklabels(species_with_component[::-1])
    
    # Set labels and title
    ax.set_xlabel('Genes', fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(f'Core iModulon {component} Comparison', fontsize=14)
    
    # Apply font to all text elements if font_path provided
    if font_path and os.path.exists(font_path):
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontproperties(font_prop)
        ax.xaxis.label.set_fontproperties(font_prop)
        ax.yaxis.label.set_fontproperties(font_prop)
        ax.title.set_fontproperties(font_prop)
    
    # Remove all spines first
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add separate boxes for left and right heatmaps
    if gap_position > 0:
        # Left heatmap box
        left_box = patches.Rectangle((0, 0), gap_position, n_species, 
                                   linewidth=1, edgecolor='black', facecolor='none')
        ax.add_patch(left_box)
        
        # Right heatmap box
        right_start = gap_position + gap_width
        right_width = total_width - right_start
        right_box = patches.Rectangle((right_start, 0), right_width, n_species,
                                    linewidth=1, edgecolor='black', facecolor='none')
        ax.add_patch(right_box)
    else:
        # Single box if no separation
        full_box = patches.Rectangle((0, 0), total_width, n_species,
                                   linewidth=1, edgecolor='black', facecolor='none')
        ax.add_patch(full_box)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#CCE5FF', edgecolor='black', linewidth=1, label='Present in iModulon'),
        Patch(facecolor='#F0F0F0', edgecolor='black', linewidth=1, label='Absent from iModulon'),
        Patch(facecolor='white', edgecolor='black', linewidth=1, label=f'Absent from {y_label.lower()}')
    ]
    
    # Position legend to the right with same gap as between heatmaps
    if gap_position > 0:
        legend_x = total_width + gap_width
    else:
        legend_x = total_width + 1
        
    # Create custom legend
    legend = ax.legend(handles=legend_elements, 
                      loc='center left',
                      bbox_to_anchor=(legend_x / total_width, 0.5),
                      frameon=True,
                      framealpha=1,
                      edgecolor='black')
    
    # Apply font to legend if provided
    if font_path and os.path.exists(font_path):
        for text in legend.get_texts():
            text.set_fontproperties(font_prop)
    
    # Adjust figure size to accommodate legend
    fig.subplots_adjust(right=0.75)
    
    # Tight layout
    plt.tight_layout()
    
    # Save or show
    if save_path:
        save_path = Path(save_path)
        if save_path.suffix in ['.svg', '.png', '.pdf', '.jpg']:
            save_file = save_path
        else:
            save_path.mkdir(parents=True, exist_ok=True)
            save_file = save_path / f"{component}_comparison_heatmap.svg"
        
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        logger.info(f"Comparison heatmap saved to {save_file}")
        
    # Always show the plot (useful for Jupyter notebooks)
    plt.show()
    
    return result_dict