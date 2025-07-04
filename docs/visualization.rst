Visualization of iModulons
==========================

This section covers the visualization functions for exploring iModulon gene weights and activities across samples.

Overview
--------

MultiModulon provides four main visualization functions:

1. **view_iModulon_weights** - Visualize gene weights within a component for a single species
2. **view_core_iModulon_weights** - Visualize a core iModulon component across all species
3. **view_iModulon_activities** - Visualize component activities across samples
4. **compare_core_iModulon_activity** - Compare core iModulon activities across multiple species for specific conditions

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

Comparing Core iModulon Activities Across Species
-------------------------------------------------

.. py:method:: MultiModulon.compare_core_iModulon_activity(component, species_in_comparison, condition_list, save_path=None, fig_size=(12, 3), font_path=None, legend_title=None, title=None)

   Compare core iModulon activities across multiple species for specific conditions.
   Creates a grouped bar plot with conditions on x-axis and species shown as different colored bars.

   :param str component: Core component name (e.g., 'Core_1', 'Core_2')
   :param list species_in_comparison: List of species names to compare
   :param list condition_list: List of conditions in format "condition:project"
   :param str save_path: Path to save the plot (optional)
   :param tuple fig_size: Figure size (default: (12, 3))
   :param str font_path: Path to custom font file (optional)
   :param str legend_title: Custom title for the legend (default: 'Species')
   :param str title: Custom title for the plot (default: 'Core iModulon {component} Activity Comparison')

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   # Compare Core_1 activities across species for specific conditions
   multiModulon.compare_core_iModulon_activity(
       component='Core_1',
       species_in_comparison=['E_coli', 'S_enterica', 'K_pneumoniae'],
       condition_list=['glucose:project1', 'lactose:project1', 'arabinose:project2']
   )

Condition Format
~~~~~~~~~~~~~~~~

Conditions must be specified as "condition:project" pairs:

.. code-block:: python

   # Comparing growth conditions from different projects
   multiModulon.compare_core_iModulon_activity(
       component='Core_1',
       species_in_comparison=['Species1', 'Species2', 'Species3'],
       condition_list=[
           'exponential:growth_study',    # Exponential phase from growth_study
           'stationary:growth_study',     # Stationary phase from growth_study
           'heat_shock:stress_project',   # Heat shock from stress_project
           'cold_shock:stress_project'    # Cold shock from stress_project
       ],
       save_path='core1_condition_comparison.svg'
   )

Understanding the Plot
~~~~~~~~~~~~~~~~~~~~~~

* **X-axis**: Conditions (grouped by the order in condition_list)
* **Y-axis**: iModulon activity values
* **Bars**: Different colors for each species
* **Dots**: Individual sample values (black dots on bars)
* **Legend**: Species names with corresponding colors

Error Handling
~~~~~~~~~~~~~~

The function validates that all conditions exist in all species:

.. code-block:: python

   # This will raise an error if any species lacks a condition
   try:
       multiModulon.compare_core_iModulon_activity(
           component='Core_1',
           species_in_comparison=['Species1', 'Species2'],
           condition_list=['rare_condition:project1']
       )
   except ValueError as e:
       print(f"Error: {e}")

Customizing Appearance
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Larger figure with custom font
   multiModulon.compare_core_iModulon_activity(
       component='Core_1',
       species_in_comparison=['Species1', 'Species2', 'Species3'],
       condition_list=['control:exp1', 'treatment:exp1'],
       fig_size=(15, 5),  # Wider figure
       font_path='/path/to/font.ttf',
       save_path='comparison_custom.svg'
   )
   
   # Custom title and legend
   multiModulon.compare_core_iModulon_activity(
       component='Core_1',
       species_in_comparison=['E_coli_K12', 'E_coli_B', 'E_coli_C'],
       condition_list=['glucose:carbon_study', 'lactose:carbon_study'],
       title='Carbon Source Response in E. coli Strains',
       legend_title='E. coli Strain',
       save_path='ecoli_carbon_response.svg'
   )

Use Cases
~~~~~~~~~

1. **Stress Response Comparison**: Compare how different species respond to the same stresses
2. **Metabolic Adaptation**: Analyze metabolic shifts across species under different carbon sources
3. **Evolutionary Analysis**: Study conservation of regulatory responses

.. code-block:: python

   # Example: Comparing stress responses
   stress_conditions = [
       'control:stress_study',
       'heat_42C:stress_study',
       'oxidative_H2O2:stress_study',
       'acid_pH5:stress_study'
   ]
   
   multiModulon.compare_core_iModulon_activity(
       component='Core_1',  # Assuming Core_1 is stress-related
       species_in_comparison=['E_coli', 'S_enterica', 'K_pneumoniae'],
       condition_list=stress_conditions,
       save_path='stress_response_comparison.svg'
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