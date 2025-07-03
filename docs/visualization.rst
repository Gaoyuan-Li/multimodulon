Visualization of iModulons
==========================

This section covers the visualization functions for exploring iModulon gene weights and activities across samples.

Overview
--------

MultiModulon provides three main visualization functions:

1. **view_iModulon_weights** - Visualize gene weights within a component for a single species
2. **view_core_iModulon_weights** - Visualize a core iModulon component across all species
3. **view_iModulon_activities** - Visualize component activities across samples

All functions support customization of appearance, highlighting, and export options.

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
   multiModulon.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       save_path='core1_weights.svg'
   )
   
   # With COG coloring
   multiModulon.view_iModulon_weights(
       species='Species1', 
       component='Core_1',
       show_COG=True,
       save_path='core1_weights_COG.svg'
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
   multiModulon.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       fig_size=(8, 6),
       font_path='/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf',
       save_path='custom_weights.svg'
   )

Visualizing Core iModulons Across Species
-----------------------------------------

.. py:method:: MultiModulon.view_core_iModulon_weights(component, save_path=None, fig_size=(6, 4), font_path=None, show_COG=False, reference_order=None)

   Visualize a core iModulon component across all species. Creates individual plots for each species
   showing the same core component, or a combined plot with subplots when COG coloring is enabled.

   :param str component: Core component name (e.g., 'Core_1', 'Core_2')
   :param str save_path: Directory path to save plots (optional)
   :param tuple fig_size: Figure size for individual plots (default: (6, 4))
   :param str font_path: Path to custom font file (optional)
   :param bool show_COG: Color genes by COG category (default: False)
   :param list reference_order: Custom species order for subplot arrangement (optional)

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   # Visualize core component across all species
   multiModulon.view_core_iModulon_weights(
       component='Core_1',
       save_path='core_plots/'
   )
   
   # With COG coloring - creates combined plot
   multiModulon.view_core_iModulon_weights(
       component='Core_1',
       show_COG=True,
       save_path='core1_all_species_COG.svg'
   )

Custom Species Order
~~~~~~~~~~~~~~~~~~~~

When using COG coloring, arrange species in a specific order:

.. code-block:: python

   # Define custom order (first 3 in top row, rest in bottom row)
   multiModulon.view_core_iModulon_weights(
       component='Core_1',
       show_COG=True,
       reference_order=['MG1655', 'BL21', 'C', 'Crooks', 'W', 'W3110'],
       save_path='core1_ordered.svg'
   )

Understanding the Output
~~~~~~~~~~~~~~~~~~~~~~~~

**Without COG coloring**: Creates individual plots for each species
   - Each plot saved as '{species}_{component}_iModulon.svg'
   - Shows gene weights on genomic coordinates
   - Includes threshold lines if available

**With COG coloring**: Creates a single combined plot
   - All species shown as subplots
   - Shared COG category legend at bottom
   - Genes colored by functional category
   - Grey dots indicate genes below threshold

Batch Processing Core Components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Plot all core components
   M = multiModulon[multiModulon.species[0]].M
   core_components = [c for c in M.columns if c.startswith('Core_')]
   
   for comp in core_components:
       # Individual species plots
       multiModulon.view_core_iModulon_weights(
           component=comp,
           save_path=f'core_plots/{comp}/'
       )
       
       # Combined COG plot
       multiModulon.view_core_iModulon_weights(
           component=comp,
           show_COG=True,
           save_path=f'core_plots/{comp}_COG.svg'
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
   multiModulon.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       save_path='core1_activities.svg'
   )
   
   # Highlight specific project
   multiModulon.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_project='ProjectA',
       save_path='core1_highlighted.svg'
   )

Condition-based Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a ``condition`` column exists in the sample sheet:

.. code-block:: python

   # Activities are averaged by condition
   # Individual sample values shown as black dots
   multiModulon.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       save_path='condition_averaged.svg'
   )
   
   # Highlight specific conditions
   multiModulon.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_condition=['Treatment1', 'Treatment2'],
       save_path='conditions_highlighted.svg'
   )

Show Only Highlighted Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Focus on specific conditions:

.. code-block:: python

   # Show only specific conditions with custom colors
   multiModulon.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_condition=['Control', 'Stress', 'Recovery'],
       show_highlight_only=True,
       show_highlight_only_color=['blue', 'red', 'green'],
       save_path='focused_conditions.svg'
   )

Multiple Highlighting Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Highlight multiple projects
   multiModulon.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_project=['ProjectA', 'ProjectB'],
       save_path='multi_project.svg'
   )
   
   # Highlight by study
   multiModulon.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_study='GSE12345',
       save_path='study_highlighted.svg'
   )

Advanced Visualization
----------------------

Batch Visualization
~~~~~~~~~~~~~~~~~~~

Create plots for multiple components:

.. code-block:: python

   # Plot all core components
   for species in multiModulon.species:
       M = multiModulon[species].M
       core_comps = [c for c in M.columns if c.startswith('Core_')]
       
       for comp in core_comps:
           # Gene weights
           multiModulon.view_iModulon_weights(
               species=species,
               component=comp,
               show_COG=True,
               save_path=f'weights/{species}_{comp}_weights.svg'
           )
           
           # Activities
           multiModulon.view_iModulon_activities(
               species=species,
               component=comp,
               save_path=f'activities/{species}_{comp}_activities.svg'
           )

Export Options
--------------

File Formats
~~~~~~~~~~~~

Save plots in different formats:

.. code-block:: python

   # Vector format (scalable)
   multiModulon.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       save_path='weights.svg'  # SVG format
   )
   
   # High-resolution raster
   multiModulon.view_iModulon_weights(
       species='Species1',
       component='Core_1',
       save_path='weights.png'  # png at 300 DPI
   )
   
   # PDF for publications
   multiModulon.view_iModulon_weights(
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
   for species in multiModulon.species:
       for comp in ['Core_1', 'Core_2', 'Unique_1']:
           # Weights without COG
           multiModulon.view_iModulon_weights(
               species=species,
               component=comp,
               save_path=f'{base_dir}/weights/{species}_{comp}.svg'
           )
           
           # Weights with COG
           multiModulon.view_iModulon_weights(
               species=species,
               component=comp,
               show_COG=True,
               save_path=f'{base_dir}/weights_COG/{species}_{comp}.svg'
           )
           
           # Activities
           multiModulon.view_iModulon_activities(
               species=species,
               component=comp,
               save_path=f'{base_dir}/activities/{species}_{comp}.svg'
           )

Best Practices
--------------

1. **Use descriptive filenames** - Include species and component names
2. **Consistent figure sizes** - Use same dimensions for comparable plots
3. **Save vector formats** - Use SVG for publication figures
4. **Document parameters** - Note thresholds and highlighting used

Next Steps
----------

1. :doc:`examples/visualization_gallery` - More visualization examples
2. Biological interpretation - Analyze visualized patterns
3. Export for further analysis - Use data in other tools