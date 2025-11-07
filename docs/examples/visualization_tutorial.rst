Characterization and Visualization 
======================

This tutorial demonstrates how to visualize and characterize iModulons after running MultiModulon analysis.

Loading Saved Results
---------------------

Start by loading a previously saved MultiModulon object:

.. code-block:: python

   from multimodulon import MultiModulon
   import pandas as pd
   import matplotlib.pyplot as plt
   
   # Load saved MultiModulon object
   multiModulon = MultiModulon.load_json_multimodulon("multiModulon_results.json.gz")

Renaming iModulons
------------------

You can rename iModulons to more meaningful names based on their function or content. There are now two separate functions for renaming core and unique iModulons:

.. code-block:: python

   # Rename a core iModulon (affects all species)
   multiModulon.rename_core_iModulon('Core_1', 'ArginineBiosynthesis')
   
   # Rename a unique iModulon (species-specific)
   multiModulon.rename_unique_iModulon('MG1655', 'Unique_1', 'StressResponse')
   
   # Note the different parameter order:
   # - rename_core_iModulon(old_name, new_name) - only 2 parameters
   # - rename_unique_iModulon(species, old_name, new_name) - 3 parameters
   
   # View the renamed components
   print(multiModulon['MG1655'].M.columns)

Important notes about renaming:

* **Core iModulons**: Use ``rename_core_iModulon()`` - automatically renames across all species where the component exists
* **Unique iModulons**: Use ``rename_unique_iModulon()`` - requires species as the first parameter
* Both methods validate that the new name doesn't already exist
* All related dataframes (M, A, M_thresholds, presence_matrix) are updated automatically
* The parameter order is different: core functions need only old and new names, while unique functions need species first

Visualizing Individual Components
---------------------------------

Gene Weights Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~

Visualize gene weights for a specific component in a species:

.. code-block:: python

   # Basic gene weight plot
   multiModulon.view_iModulon_weights(
       species="Species1",
       component="Core_1",
       fig_size=(6, 4),
       font_path="/path/to/font.ttf"  # optional
   )
   
   # With COG category coloring
   multiModulon.view_iModulon_weights(
       species="Species1",
       component="Core_1",
       fig_size=(6, 4),
       show_COG=True
   )

The plot shows:

* Genes sorted by weight magnitude
* Red threshold line (if optimized)
* COG categories in color (if enabled)
* Top gene names on the right

Viewing Component Genes
~~~~~~~~~~~~~~~~~~~~~~~

Get detailed information about genes in a component:

.. code-block:: python

   # View gene table for a component
   gene_info = multiModulon.view_iModulon_genes(
       species="Species1",
       component="Core_1"
   )

This returns a subset of the gene table.

Activity Visualization
~~~~~~~~~~~~~~~~~~~~~~

Visualize component activities across samples:

.. code-block:: python

   # Basic activity plot
   multiModulon.view_iModulon_activities(
       species="Species1",
       component="Core_1",
       fig_size=(12, 3)
   )
   
   # Highlight specific projects
   multiModulon.view_iModulon_activities(
       species="Species1",
       component="Core_1",
       highlight_project="ProjectA"
   )
   
   # Highlight multiple projects
   multiModulon.view_iModulon_activities(
       species="Species1",
       component="Core_1",
       highlight_project=["ProjectA", "ProjectB"]
   )

Features:

* Bar plot of activities
* Project/study grouping on x-axis
* Color highlighting for specific projects

Comparing Core Components Across Species
----------------------------------------

Visualize how core components are conserved across species:

.. code-block:: python

   # Compare gene weights across species
   multiModulon.view_core_iModulon_weights(
       component="Core_1",
       fig_size=(6, 4),
       reference_order=['Species1', 'Species2', 'Species3'],
       show_COG=True
   )

Gene Membership Comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~

Create detailed comparison of gene membership across species:

.. code-block:: python

   # Generate membership comparison
   comparison_df = multiModulon.compare_core_iModulon(
       component='Core_1',
       y_label='Strains',
       reference_order=['Species1', 'Species2', 'Species3'],
       fig_size=(20, 6),
       font_path="/path/to/font.ttf",
       save_path="output_dir/"
   )

This creates:

* Heatmap showing gene presence across species
* Genes grouped by conservation pattern
* Visual identification of core vs species-specific genes

To highlight only a subset of genes on the heatmap, pass the new ``show_list`` parameter
with x-tick labels (locus tags by default, gene names when ``show_gene_names=True``):

.. code-block:: python

   multiModulon.compare_core_iModulon(
       component='Core_1',
       show_gene_names=True,
       show_list=['cysB', 'cysK', 'cysM']
   )

Characterizing Unique Components
--------------------------------

Explore species-specific regulatory modules:

.. code-block:: python

   # Visualize unique component for a species
   multiModulon.view_iModulon_weights(
       species="Species1",
       component="Unique_1",
       fig_size=(6, 4),
       show_COG=True
   )
   
   # Check activities
   multiModulon.view_iModulon_activities(
       species="Species1",
       component="Unique_1",
       fig_size=(12, 3)
   )

Advanced Visualization Options
------------------------------

Condition-based Analysis
~~~~~~~~~~~~~~~~~~~~~~~~

When sample sheet contains a "condition" column:

.. code-block:: python

   # Activities are automatically averaged by condition
   multiModulon.view_iModulon_activities(
       species="Species1",
       component="Core_1",
       highlight_condition=["Control", "Treatment"]
   )
   
   # Show only specific conditions
   multiModulon.view_iModulon_activities(
       species="Species1",
       component="Core_1",
       highlight_condition=["Control", "Stress", "Recovery"],
       show_highlight_only=True,
       show_highlight_only_color=["blue", "red", "green"]
   )

Custom Styling
~~~~~~~~~~~~~~

Customize plot appearance:

.. code-block:: python

   # Custom figure size and font
   multiModulon.view_iModulon_weights(
       species="Species1",
       component="Core_1",
       fig_size=(8, 6),
       font_path="/usr/share/fonts/truetype/arial.ttf",
       save_path="custom_plot.svg"  
   )

Batch Visualization
~~~~~~~~~~~~~~~~~~~

Visualize all components systematically:

.. code-block:: python

   # Get all components for a species
   M = multiModulon['Species1'].M
   components = M.columns
   
   # Separate core and unique
   core_components = [c for c in components if c.startswith('Core_')]
   unique_components = [c for c in components if c.startswith('Unique_')]
   
   # Batch visualize
   for comp in core_components:
       multiModulon.view_iModulon_weights(
           species="Species1",
           component=comp,
           show_COG=True,
           save_path=f"weights/{comp}_weights.svg"
       )
       
       multiModulon.view_iModulon_activities(
           species="Species1",
           component=comp,
           save_path=f"activities/{comp}_activities.svg"
       )

Interpreting Results
--------------------

Core Components
~~~~~~~~~~~~~~~

Core components represent conserved regulatory modules:

* High conservation across species indicates fundamental regulation
* Differences in gene membership reveal species adaptations

Unique Components
~~~~~~~~~~~~~~~~~

Unique components capture species-specific regulation:

* May represent adaptation to specific environments
* Could indicate gain/loss of regulatory mechanisms

Export for Further Analysis
---------------------------

Export data for external tools:

.. code-block:: python

   # Export component genes
   for comp in core_components:
       genes = multiModulon.view_iModulon_genes("Species1", comp)
       genes.to_csv(f"{comp}_genes.csv")
   
   # Export activities
   A = multiModulon['Species1'].A
   A.to_csv("Species1_activities.csv")
   
   # Export for gene set enrichment
   presence = multiModulon['Species1'].presence_matrix
   for comp in presence.columns:
       gene_list = presence[presence[comp] == 1].index
       with open(f"{comp}_genelist.txt", 'w') as f:
           f.write('\n'.join(gene_list))
