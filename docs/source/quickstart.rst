Quick Start Guide
=================

This guide will walk you through a basic MultiModulon analysis workflow.

Overview
--------

A typical MultiModulon analysis consists of these steps:

1. Prepare your data in the required format
2. Initialize the MultiModulon object
3. Generate BBH and align genes across species
4. Run multi-view ICA to identify regulatory modules
5. Visualize and interpret results

Basic Workflow
--------------

Here's a complete example workflow:

.. code-block:: python

   from multimodulon import MultiModulon
   
   # Step 1: Initialize MultiModulon with your data directory
   mm = MultiModulon("path/to/Input_Data")
   
   # View loaded species
   mm.summary()
   
   # Step 2: Generate Bidirectional Best Hits (BBH)
   mm.generate_BBH(output_path="Output_BBH", threads=8)
   
   # Step 3: Align genes across species
   combined_gene_db = mm.align_genes(
       input_bbh_dir="Output_BBH",
       output_dir="Output_Gene_Info",
       bbh_threshold=40  # minimum percent identity
   )
   
   # Step 4: Generate expression matrices
   mm.generate_X("Output_Gene_Info")
   
   # Step 5: Optimize component numbers (optional but recommended)
   optimal_core, scores = mm.optimize_number_of_core_components(
       max_k=30,
       step=5,
       save_plot="optimization_plot.png"
   )
   
   optimal_unique, optimal_total = mm.optimize_number_of_unique_components(
       optimal_num_core_components=optimal_core,
       step=5,
       save_plots="optimization_plots/"
   )
   
   # Step 6: Run multi-view ICA
   M_matrices, A_matrices = mm.run_robust_multiview_ica(
       a=optimal_total,  # e.g., {'species1': 50, 'species2': 60}
       c=optimal_core,   # e.g., 20
       num_runs=100,
       mode='gpu'
   )
   
   # Step 7: Generate activity matrices
   mm.generate_A()
   
   # Step 8: Visualize results
   # View iModulon gene weights
   mm.view_iModulon_weights(
       species='species1',
       component='Core_1',
       save_path='core1_weights.png',
       show_COG=True
   )
   
   # View iModulon activities across samples
   mm.view_iModulon_activities(
       species='species1',
       component='Core_1',
       save_path='core1_activities.png',
       highlight_project='ProjectA'
   )

Working with Results
--------------------

Access specific data:

.. code-block:: python

   # Get data for a specific species
   species_data = mm['species1']
   
   # Access matrices
   X_matrix = species_data.X  # Expression matrix
   M_matrix = species_data.M  # Mixing matrix (gene weights)
   A_matrix = species_data.A  # Activity matrix
   
   # Access metadata
   sample_sheet = species_data.sample_sheet
   gene_table = species_data.gene_table

Export results:

.. code-block:: python

   # Save M matrices
   for species in mm.species:
       M = mm[species].M
       M.to_csv(f"{species}_M_matrix.csv")
   
   # Save A matrices
   for species in mm.species:
       A = mm[species].A
       A.to_csv(f"{species}_A_matrix.csv")

Next Steps
----------

* See :doc:`data_preparation` for detailed data format requirements
* Check :doc:`optimization` for advanced optimization strategies
* Explore :doc:`visualization` for all visualization options
* Read :doc:`examples/basic_workflow` for a complete worked example