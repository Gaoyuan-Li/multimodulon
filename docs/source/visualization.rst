Visualization of iModulons
==========================

This section covers the visualization functions for exploring iModulon gene weights and activities across samples.

Overview
--------

MultiModulon provides two main visualization functions:

1. **view_iModulon_weights** - Visualize gene weights within a component
2. **view_iModulon_activities** - Visualize component activities across samples

Both functions support customization of appearance, highlighting, and export options.

Visualizing Gene Weights
------------------------

.. py:method:: MultiModulon.view_iModulon_weights(species, component, save_path=None, fig_size=(6, 4), font_path=None, show_COG=False)

   Create a bar plot showing gene weights for a specific iModulon component.

   :param str species: Species/strain name
   :param str component: Component name (e.g., 'Core_1', 'Unique_1')
   :param str save_path: Path to save the plot (optional)
   :param tuple fig_size: Figure size as (width, height) (default: (6, 4))
   :param str font_path: Path to custom font file (optional)
   :param bool show_COG: Color genes by COG category (default: False)

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   # Simple gene weight plot
   mm.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       save_path='core1_weights.png'
   )
   
   # With COG coloring
   mm.view_iModulon_weights(
       species='Species1', 
       component='Core_1',
       show_COG=True,
       save_path='core1_weights_COG.png'
   )

Understanding the Plot
~~~~~~~~~~~~~~~~~~~~~~

* **X-axis**: Genes sorted by weight magnitude
* **Y-axis**: Gene weights (coefficients from M matrix)
* **Red line**: Threshold (if optimized)
* **Colors**: COG categories (if show_COG=True)
* **Labels**: Top genes shown on right side

COG Categories
~~~~~~~~~~~~~~

When ``show_COG=True``, genes are colored by functional category:

.. code-block:: python

   # COG categories and their colors:
   # - Translation (J): black
   # - Transcription (K): sandybrown  
   # - Replication (L): fuchsia
   # - Cell division (D): olive
   # - Defense (V): orchid
   # - Signal transduction (T): teal
   # - Cell membrane (M): purple
   # - Energy production (C): red
   # - Carbohydrate metabolism (G): gold
   # - Amino acid metabolism (E): darkgreen
   # - Nucleotide metabolism (F): pink
   # - Coenzyme metabolism (H): brown
   # - Lipid metabolism (I): lightsalmon
   # - Inorganic ion metabolism (P): darkblue
   # - Secondary metabolism (Q): sienna
   # - Unknown function (S): lightgray
   # - Not in COG: gray

Customizing Appearance
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Larger figure with custom font
   mm.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       fig_size=(8, 6),
       font_path='/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf',
       save_path='custom_weights.png'
   )

Visualizing iModulon Activities
-------------------------------

.. py:method:: MultiModulon.view_iModulon_activities(species, component, save_path=None, fig_size=(12, 3), font_path=None, highlight_project=None, highlight_study=None, highlight_condition=None, show_highlight_only=False, show_highlight_only_color=None)

   Create a bar plot showing component activities across samples.

   :param str species: Species/strain name
   :param str component: Component name
   :param str save_path: Path to save the plot
   :param tuple fig_size: Figure size (default: (12, 3))
   :param str font_path: Path to custom font
   :param highlight_project: Project(s) to highlight (str or list)
   :param str highlight_study: Study to highlight
   :param highlight_condition: Condition(s) to highlight (str or list)
   :param bool show_highlight_only: Only show highlighted conditions
   :param list show_highlight_only_color: Colors for highlighted conditions

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   # Simple activity plot
   mm.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       save_path='core1_activities.png'
   )
   
   # Highlight specific project
   mm.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_project='ProjectA',
       save_path='core1_highlighted.png'
   )

Condition-based Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a ``condition`` column exists in the sample sheet:

.. code-block:: python

   # Activities are averaged by condition
   # Individual sample values shown as black dots
   mm.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       save_path='condition_averaged.png'
   )
   
   # Highlight specific conditions
   mm.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_condition=['Treatment1', 'Treatment2'],
       save_path='conditions_highlighted.png'
   )

Show Only Highlighted Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Focus on specific conditions:

.. code-block:: python

   # Show only specific conditions with custom colors
   mm.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_condition=['Control', 'Stress', 'Recovery'],
       show_highlight_only=True,
       show_highlight_only_color=['blue', 'red', 'green'],
       save_path='focused_conditions.png'
   )

Multiple Highlighting Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Highlight multiple projects
   mm.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_project=['ProjectA', 'ProjectB'],
       save_path='multi_project.png'
   )
   
   # Highlight by study
   mm.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_study='GSE12345',
       save_path='study_highlighted.png'
   )

Advanced Visualization
----------------------

Batch Visualization
~~~~~~~~~~~~~~~~~~~

Create plots for multiple components:

.. code-block:: python

   # Plot all core components
   for species in mm.species:
       M = mm[species].M
       core_comps = [c for c in M.columns if c.startswith('Core_')]
       
       for comp in core_comps:
           # Gene weights
           mm.view_iModulon_weights(
               species=species,
               component=comp,
               show_COG=True,
               save_path=f'weights/{species}_{comp}_weights.png'
           )
           
           # Activities
           mm.view_iModulon_activities(
               species=species,
               component=comp,
               save_path=f'activities/{species}_{comp}_activities.png'
           )

Custom Plotting
~~~~~~~~~~~~~~~

Access the data directly for custom plots:

.. code-block:: python

   import matplotlib.pyplot as plt
   import seaborn as sns
   
   # Get component data
   species = 'Species1'
   component = 'Core_1'
   
   # Gene weights
   M = mm[species].M
   weights = M[component].sort_values(ascending=False)
   
   # Custom weight plot
   plt.figure(figsize=(10, 6))
   plt.bar(range(len(weights)), weights.values)
   plt.xlabel('Genes (sorted by weight)')
   plt.ylabel('Weight')
   plt.title(f'{component} gene weights in {species}')
   
   # Annotate top genes
   for i, (gene, weight) in enumerate(weights.head(5).items()):
       plt.annotate(gene, (i, weight), rotation=45)
   
   plt.tight_layout()
   plt.savefig('custom_weights.png')

Activity Heatmaps
~~~~~~~~~~~~~~~~~

Create heatmaps for multiple components:

.. code-block:: python

   # Get activity matrix
   A = mm['Species1'].A
   
   # Select components of interest
   components = ['Core_1', 'Core_2', 'Core_3', 'Unique_1', 'Unique_2']
   A_subset = A.loc[components]
   
   # Create heatmap
   plt.figure(figsize=(15, 5))
   sns.heatmap(
       A_subset,
       cmap='RdBu_r',
       center=0,
       cbar_kws={'label': 'Activity'},
       yticklabels=True
   )
   plt.xlabel('Samples')
   plt.ylabel('Components')
   plt.title('iModulon activity heatmap')
   plt.tight_layout()
   plt.savefig('activity_heatmap.png')

Comparative Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~

Compare components across species:

.. code-block:: python

   # Compare Core_1 across species
   fig, axes = plt.subplots(len(mm.species), 1, figsize=(12, 4*len(mm.species)))
   
   component = 'Core_1'
   for i, species in enumerate(mm.species):
       ax = axes[i] if len(mm.species) > 1 else axes
       
       # Get activities
       A = mm[species].A
       activities = A.loc[component]
       
       # Plot
       ax.bar(range(len(activities)), activities.values)
       ax.set_title(f'{component} activities in {species}')
       ax.set_xlabel('Samples')
       ax.set_ylabel('Activity')
       ax.axhline(y=0, color='black', linewidth=0.5)
   
   plt.tight_layout()
   plt.savefig('core1_comparison.png')

Export Options
--------------

File Formats
~~~~~~~~~~~~

Save plots in different formats:

.. code-block:: python

   # Vector format (scalable)
   mm.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       save_path='weights.svg'  # SVG format
   )
   
   # High-resolution raster
   mm.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       save_path='weights.png'  # PNG at 300 DPI
   )
   
   # PDF for publications
   mm.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       save_path='weights.pdf'
   )

Directory Organization
~~~~~~~~~~~~~~~~~~~~~~

Organize outputs systematically:

.. code-block:: python

   import os
   
   # Create directory structure
   base_dir = 'imodulon_plots'
   for subdir in ['weights', 'activities', 'weights_COG']:
       os.makedirs(f'{base_dir}/{subdir}', exist_ok=True)
   
   # Save with organized naming
   for species in mm.species:
       for comp in ['Core_1', 'Core_2', 'Unique_1']:
           # Weights without COG
           mm.view_iModulon_weights(
               species=species,
               component=comp,
               save_path=f'{base_dir}/weights/{species}_{comp}.png'
           )
           
           # Weights with COG
           mm.view_iModulon_weights(
               species=species,
               component=comp,
               show_COG=True,
               save_path=f'{base_dir}/weights_COG/{species}_{comp}.png'
           )
           
           # Activities
           mm.view_iModulon_activities(
               species=species,
               component=comp,
               save_path=f'{base_dir}/activities/{species}_{comp}.png'
           )

Best Practices
--------------

1. **Use descriptive filenames** - Include species and component names
2. **Consistent figure sizes** - Use same dimensions for comparable plots
3. **Save vector formats** - Use SVG/PDF for publication figures
4. **Document parameters** - Note thresholds and highlighting used
5. **Check orientations** - Ensure text is readable in saved files

Troubleshooting
---------------

**Font warnings:**

.. code-block:: python

   # Use system fonts or specify path
   import matplotlib.font_manager as fm
   
   # List available fonts
   fonts = fm.findSystemFonts()
   print("Available fonts:", fonts[:5])
   
   # Use a specific font
   mm.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       font_path=fonts[0]  # Use first available
   )

**Large datasets (many samples):**

.. code-block:: python

   # Increase figure width for many samples
   n_samples = len(mm['Species1'].A.columns)
   fig_width = max(12, n_samples * 0.1)  # Scale with samples
   
   mm.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       fig_size=(fig_width, 3)
   )

**Memory issues with batch plotting:**

.. code-block:: python

   # Close figures after saving
   import matplotlib.pyplot as plt
   
   for comp in components:
       mm.view_iModulon_weights(
           species='Species1',
           component=comp,
           save_path=f'{comp}.png'
       )
       plt.close('all')  # Free memory

Next Steps
----------

1. :doc:`examples/visualization_gallery` - More visualization examples
2. Biological interpretation - Analyze visualized patterns
3. Export for further analysis - Use data in other tools