Visualization API Reference
===========================

This section provides detailed API documentation for visualization functions.

Gene Weight Visualization
-------------------------

.. automethod:: multimodulon.MultiModulon.view_iModulon_weights

   Visualize gene weights for a specific iModulon component.

   **Parameters:**
   
   * **species** (*str*) -- Species/strain name
   * **component** (*str*) -- Component name (e.g., 'Core_1', 'Unique_1')
   * **save_path** (*str, optional*) -- Path to save plot. Can be:
     
     - Full file path: ``'/path/to/plot.png'``
     - Directory path: ``'/path/to/dir/'`` (saves as ``{species}_{component}_weights.svg``)
     - ``None``: Display plot without saving
   
   * **fig_size** (*tuple*) -- Figure size as (width, height) (default: (6, 4))
   * **font_path** (*str, optional*) -- Path to font file (e.g., Arial.ttf)
   * **show_COG** (*bool*) -- Color genes by COG category (default: False)
   
   **Raises:**
   
   * **ValueError** -- If species not found, M matrix not available, or component not found
   
   **Plot Elements:**
   
   * Bar plot of gene weights sorted by magnitude
   * Red dashed line indicating threshold (if available)
   * Gene names for top weighted genes
   * COG category coloring (if enabled)
   * Legend showing COG categories (if enabled)
   
   **Example:**
   
   .. code-block:: python
      
      # Basic usage
      mm.view_iModulon_weights(
          species='E_coli',
          component='Core_1',
          save_path='core1_weights.png'
      )
      
      # With COG coloring and custom appearance
      mm.view_iModulon_weights(
          species='E_coli',
          component='Core_1',
          show_COG=True,
          fig_size=(8, 6),
          font_path='/usr/share/fonts/truetype/arial.ttf',
          save_path='core1_weights_COG.svg'
      )

Activity Visualization  
----------------------

.. automethod:: multimodulon.MultiModulon.view_iModulon_activities

   Visualize iModulon activities across samples with advanced grouping and highlighting.

   **Parameters:**
   
   * **species** (*str*) -- Species/strain name
   * **component** (*str*) -- Component name
   * **save_path** (*str, optional*) -- Path to save plot
   * **fig_size** (*tuple*) -- Figure size (default: (12, 3))
   * **font_path** (*str, optional*) -- Path to font file
   * **highlight_project** (*str or list, optional*) -- Project(s) to highlight
   * **highlight_study** (*str, optional*) -- Study accession to highlight
   * **highlight_condition** (*str or list, optional*) -- Condition(s) to highlight
   * **show_highlight_only** (*bool*) -- Only show highlighted conditions (default: False)
   * **show_highlight_only_color** (*list, optional*) -- Colors for highlighted conditions
   
   **Behavior Modes:**
   
   1. **Default mode**: Shows all samples as individual bars
   2. **Condition mode**: When 'condition' column exists, shows averaged bars with dots
   3. **Highlight mode**: Colors specific projects/studies/conditions
   4. **Focus mode**: When ``show_highlight_only=True``, shows only selected conditions
   
   **Example:**
   
   .. code-block:: python
      
      # Basic activity plot
      mm.view_iModulon_activities(
          species='E_coli',
          component='Core_1'
      )
      
      # Highlight specific projects
      mm.view_iModulon_activities(
          species='E_coli',
          component='Core_1',
          highlight_project=['ProjectA', 'ProjectB'],
          save_path='highlighted_activities.png'
      )
      
      # Focus on specific conditions
      mm.view_iModulon_activities(
          species='E_coli',
          component='Core_1',
          highlight_condition=['Control', 'Stress', 'Recovery'],
          show_highlight_only=True,
          show_highlight_only_color=['blue', 'red', 'green'],
          save_path='focused_conditions.png'
      )

COG Categories
--------------

The following COG categories are used for coloring genes:

.. list-table:: COG Categories and Colors
   :widths: 40 40 20
   :header-rows: 1
   
   * - Category
     - Description
     - Color
   * - J
     - Translation, ribosomal structure
     - black
   * - K
     - Transcription
     - sandybrown
   * - L
     - Replication, recombination, repair
     - fuchsia
   * - D
     - Cell cycle control, division
     - olive
   * - V
     - Defense mechanisms
     - orchid
   * - T
     - Signal transduction
     - teal
   * - M
     - Cell wall/membrane biogenesis
     - purple
   * - N
     - Cell motility
     - orange
   * - U
     - Intracellular trafficking
     - cyan
   * - O
     - Protein turnover, chaperones
     - yellow
   * - X
     - Mobilome
     - lime
   * - C
     - Energy production
     - red
   * - G
     - Carbohydrate metabolism
     - gold
   * - E
     - Amino acid metabolism
     - darkgreen
   * - F
     - Nucleotide metabolism
     - pink
   * - H
     - Coenzyme metabolism
     - brown
   * - I
     - Lipid metabolism
     - lightsalmon
   * - P
     - Inorganic ion metabolism
     - darkblue
   * - Q
     - Secondary metabolites
     - sienna
   * - R
     - General function
     - darkgray
   * - S
     - Unknown function
     - lightgray
   * - (not in COG)
     - No COG assignment
     - gray

Advanced Visualization Examples
-------------------------------

Creating Publication Figures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import matplotlib.pyplot as plt
   from matplotlib import rcParams
   
   # Set publication parameters
   rcParams['font.family'] = 'sans-serif'
   rcParams['font.sans-serif'] = ['Arial']
   rcParams['font.size'] = 12
   
   # Create high-quality figures
   for component in ['Core_1', 'Core_2', 'Core_3']:
       # Gene weights with COG
       mm.view_iModulon_weights(
           species='E_coli',
           component=component,
           show_COG=True,
           fig_size=(8, 6),
           save_path=f'figures/{component}_weights.pdf'
       )
       
       # Activities with conditions
       mm.view_iModulon_activities(
           species='E_coli', 
           component=component,
           fig_size=(14, 4),
           save_path=f'figures/{component}_activities.pdf'
       )

Batch Visualization
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import os
   
   # Create output directories
   os.makedirs('plots/weights', exist_ok=True)
   os.makedirs('plots/activities', exist_ok=True)
   
   # Visualize all components
   for species in mm.species:
       M = mm[species].M
       if M is None:
           continue
           
       for component in M.columns:
           # Skip if below threshold
           if M[component].abs().max() < 2:
               continue
           
           # Gene weights
           try:
               mm.view_iModulon_weights(
                   species=species,
                   component=component,
                   show_COG=True,
                   save_path=f'plots/weights/{species}_{component}.png'
               )
           except Exception as e:
               print(f"Error plotting {species} {component}: {e}")
           
           # Activities
           try:
               mm.view_iModulon_activities(
                   species=species,
                   component=component,
                   save_path=f'plots/activities/{species}_{component}.png'
               )
           except Exception as e:
               print(f"Error plotting {species} {component}: {e}")

Custom Styling
~~~~~~~~~~~~~~

.. code-block:: python

   # Modify plot after creation
   import matplotlib.pyplot as plt
   
   # Create plot but don't save
   mm.view_iModulon_weights(
       species='E_coli',
       component='Core_1',
       show_COG=True
   )
   
   # Get current figure and axis
   fig = plt.gcf()
   ax = plt.gca()
   
   # Customize
   ax.set_title('Core Module 1: Ribosomal Genes', fontsize=16, fontweight='bold')
   ax.set_xlabel('Genes sorted by weight', fontsize=14)
   ax.set_ylabel('Gene weight', fontsize=14)
   
   # Add annotations
   ax.annotate('Ribosomal proteins', 
               xy=(10, 15), xytext=(20, 20),
               arrowprops=dict(arrowstyle='->', color='red'),
               fontsize=12, color='red')
   
   # Save with tight layout
   plt.tight_layout()
   plt.savefig('custom_core1.pdf', dpi=300, bbox_inches='tight')

Comparative Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Compare same component across species
   fig, axes = plt.subplots(2, 2, figsize=(12, 10))
   axes = axes.flatten()
   
   component = 'Core_1'
   species_list = mm.species[:4]  # First 4 species
   
   for i, species in enumerate(species_list):
       # Get M matrix and component
       M = mm[species].M
       if M is None or component not in M.columns:
           axes[i].text(0.5, 0.5, 'No data', ha='center', va='center')
           axes[i].set_title(species)
           continue
       
       # Plot gene weights
       weights = M[component].sort_values(ascending=False)
       axes[i].bar(range(len(weights)), weights.values, color='lightblue')
       axes[i].axhline(y=0, color='black', linewidth=0.5)
       axes[i].set_title(f'{species}: {component}')
       axes[i].set_xlabel('Genes')
       axes[i].set_ylabel('Weight')
       
       # Add threshold if available
       if hasattr(mm[species], 'M_thresholds') and mm[species].M_thresholds:
           threshold = mm[species].M_thresholds.get(component, 0)
           axes[i].axhline(y=threshold, color='red', linestyle='--', alpha=0.7)
           axes[i].axhline(y=-threshold, color='red', linestyle='--', alpha=0.7)
   
   plt.tight_layout()
   plt.savefig('core1_comparison.png', dpi=300)

Interactive Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Create interactive plot with plotly
   import plotly.graph_objects as go
   from plotly.subplots import make_subplots
   
   # Get data
   species = 'E_coli'
   component = 'Core_1'
   M = mm[species].M
   A = mm[species].A
   
   # Create subplots
   fig = make_subplots(
       rows=2, cols=1,
       subplot_titles=('Gene Weights', 'Sample Activities'),
       vertical_spacing=0.12
   )
   
   # Gene weights
   weights = M[component].sort_values(ascending=False)
   fig.add_trace(
       go.Bar(x=list(range(len(weights))), 
              y=weights.values,
              text=weights.index,
              hovertemplate='Gene: %{text}<br>Weight: %{y:.2f}'),
       row=1, col=1
   )
   
   # Activities
   activities = A.loc[component]
   fig.add_trace(
       go.Bar(x=list(range(len(activities))),
              y=activities.values,
              text=activities.index,
              hovertemplate='Sample: %{text}<br>Activity: %{y:.2f}'),
       row=2, col=1
   )
   
   # Update layout
   fig.update_layout(
       title=f'{species} - {component}',
       showlegend=False,
       height=800
   )
   
   # Save interactive HTML
   fig.write_html(f'{species}_{component}_interactive.html')