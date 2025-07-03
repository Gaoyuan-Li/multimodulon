Visualization Tutorial
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

This returns a DataFrame with:

* Gene identifiers
* Gene names and products
* COG categories
* Gene weights
* Presence in component (0/1)

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
* Condition averaging (if applicable)
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
       save_path="custom_plot.pdf"
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
           save_path=f"weights/{comp}_weights.png"
       )
       
       multiModulon.view_iModulon_activities(
           species="Species1",
           component=comp,
           save_path=f"activities/{comp}_activities.png"
       )

Interpreting Results
--------------------

Core Components
~~~~~~~~~~~~~~~

Core components represent conserved regulatory modules:

* High conservation across species indicates fundamental regulation
* Differences in gene membership reveal species adaptations
* Activity patterns show condition-specific regulation

Unique Components
~~~~~~~~~~~~~~~~~

Unique components capture species-specific regulation:

* May represent adaptation to specific environments
* Could indicate gain/loss of regulatory mechanisms
* Activities often correlate with species-specific conditions

Quality Assessment
~~~~~~~~~~~~~~~~~~

Check component quality:

.. code-block:: python

   # Access thresholds
   thresholds = multiModulon['Species1'].M_thresholds
   
   # Check genes per component
   presence = multiModulon['Species1'].presence_matrix
   genes_per_component = presence.sum()
   
   print("Genes per component:")
   print(genes_per_component.sort_values(ascending=False))

Best Practices
--------------

1. **Start with core components** - These are most reliable
2. **Check gene enrichment** - Look for functional coherence
3. **Examine activities** - Should correlate with known biology
4. **Compare across species** - Conservation validates findings
5. **Document interpretations** - Save plots with annotations

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