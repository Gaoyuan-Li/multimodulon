Core API Reference
==================

This section provides detailed API documentation for the core MultiModulon classes.

MultiModulon Class
------------------

.. autoclass:: multimodulon.MultiModulon
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __getitem__
   
   .. rubric:: Core Methods
   
   .. automethod:: __init__
   .. automethod:: __getitem__
   .. automethod:: summary
   
   .. rubric:: Properties
   
   .. autoproperty:: species
   .. autoproperty:: bbh
   .. autoproperty:: combined_gene_db
   
   .. rubric:: BBH and Gene Alignment
   
   .. automethod:: generate_BBH
   .. automethod:: align_genes
   .. automethod:: save_bbh
   .. automethod:: load_bbh
   .. automethod:: get_orthologs
   
   .. rubric:: Expression Matrix Methods
   
   .. automethod:: generate_X
   .. automethod:: generate_A
   
   .. rubric:: ICA Methods
   
   .. automethod:: run_multiview_ica
   .. automethod:: run_robust_multiview_ica
   
   .. rubric:: Optimization Methods
   
   .. automethod:: optimize_number_of_core_components
   .. automethod:: optimize_number_of_unique_components
   
   .. rubric:: Analysis Methods
   
   .. automethod:: optimize_M_thresholds
   .. automethod:: calculate_explained_variance
   
   .. rubric:: Visualization Methods
   
   .. automethod:: view_iModulon_weights
   .. automethod:: view_iModulon_activities
   
   .. rubric:: Utility Methods
   
   .. automethod:: gff2pandas
   .. automethod:: create_gene_table
   .. automethod:: add_eggnog_annotation

SpeciesData Class
-----------------

.. autoclass:: multimodulon.species_data.SpeciesData
   :members:
   :undoc-members:
   :show-inheritance:
   
   .. rubric:: Properties (Lazy-loaded)
   
   .. autoproperty:: log_tpm
   .. autoproperty:: log_tpm_norm
   .. autoproperty:: X
   .. autoproperty:: M
   .. autoproperty:: A
   .. autoproperty:: sample_sheet
   .. autoproperty:: gene_table
   .. autoproperty:: M_thresholds
   .. autoproperty:: presence_matrix
   
   .. rubric:: Methods
   
   .. automethod:: validate_data

Example Usage
-------------

Basic Workflow
~~~~~~~~~~~~~~

.. code-block:: python

   from multimodulon import MultiModulon
   
   # Initialize
   mm = MultiModulon("path/to/Input_Data")
   
   # Generate BBH and align genes
   mm.generate_BBH(threads=8)
   combined_gene_db = mm.align_genes(bbh_threshold=40)
   
   # Generate expression matrices
   mm.generate_X("Output_Gene_Info")
   
   # Run optimization
   optimal_core, _ = mm.optimize_number_of_core_components(
       max_k=30,
       step=5,
       metric='effect_size'
   )
   
   optimal_unique, optimal_total = mm.optimize_number_of_unique_components(
       optimal_num_core_components=optimal_core,
       step=5
   )
   
   # Run multi-view ICA
   M_matrices, A_matrices = mm.run_robust_multiview_ica(
       a=optimal_total,
       c=optimal_core,
       num_runs=100
   )
   
   # Visualize results
   mm.view_iModulon_weights('Species1', 'Core_1', show_COG=True)
   mm.view_iModulon_activities('Species1', 'Core_1')

Accessing Data
~~~~~~~~~~~~~~

.. code-block:: python

   # Access species data
   species_data = mm['Species1']
   
   # Get matrices
   expression = species_data.log_tpm
   gene_weights = species_data.M
   activities = species_data.A
   
   # Get metadata
   samples = species_data.sample_sheet
   genes = species_data.gene_table
   
   # Check data availability
   if species_data.M is not None:
       print("ICA has been performed")
   else:
       print("ICA not yet performed")

Advanced Options
~~~~~~~~~~~~~~~~

.. code-block:: python

   # Custom BBH analysis
   from multimodulon.utils.bbh import BBHAnalyzer
   
   analyzer = BBHAnalyzer(
       docker_image="custom/blast:latest"
   )
   
   # Manual threshold optimization
   species_data = mm['Species1']
   M = species_data.M
   
   # Apply custom threshold
   threshold = 2.5
   presence = (M.abs() > threshold).astype(int)
   
   # Export results
   M.to_csv("Species1_M_matrix.csv")
   A.to_csv("Species1_A_matrix.csv")