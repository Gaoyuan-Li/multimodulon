Basic Workflow
==============

This example demonstrates the complete MultiModulon analysis workflow from data loading to component visualization.

Step 1: Initialize MultiModulon Object
--------------------------------------

Load data from the Input_Data directory containing expression matrices, gene annotations, and sample metadata for all strains.

.. code-block:: python

   from multimodulon import MultiModulon
   import pandas as pd
   import numpy as np
   
   # Initialize MultiModulon object
   multiModulon = MultiModulon(input_data_path)

The initialization will:

* Scan the input directory for species/strain subdirectories
* Load expression data (log_tpm.csv)
* Load sample metadata (sample_sheet.csv)
* Validate data consistency
* Report the number of genes and samples for each species

Step 2: Generate BBH Files
--------------------------

Generate Bidirectional Best Hits (BBH) files for ortholog detection between all strain pairs.

.. code-block:: python

   # Generate BBH files using multiple threads
   multiModulon.generate_BBH(output_path, threads=8)

This step:

* Extracts protein sequences from genome files
* Runs BLAST between all species pairs
* Identifies reciprocal best hits
* Saves results as CSV files

Step 3: Align Genes Across Strains
-----------------------------------

Create a unified gene database by aligning genes across all strains using the BBH results.

.. code-block:: python

   # Align genes across all strains
   combined_gene_db = multiModulon.align_genes(
       input_bbh_dir=output_bbh_path,
       output_dir=output_gene_info_path,
       reference_order=['Species1', 'Species2', 'Species3'],  # optional
       bbh_threshold=90  # optional: minimum percent identity
   )

The alignment process:

* Groups orthologous genes using Union-Find algorithm
* Creates gene families across species
* Generates a combined gene database
* Saves the alignment for future use

Step 4: Create Gene Tables
--------------------------

Parse GFF files to create gene annotation tables for each strain.

.. code-block:: python

   # Create gene tables from GFF files
   multiModulon.create_gene_table()
   
   # Optional: Add eggNOG annotations
   multiModulon.add_eggnog_annotation(eggnog_output_path)

This enriches gene information with:

* Gene names and products
* COG categories
* Functional annotations

Step 5: Generate Aligned Expression Matrices
--------------------------------------------

Create expression matrices with consistent gene indexing across all strains.

.. code-block:: python

   # Generate aligned expression matrices
   multiModulon.generate_X(gene_info_folder_path)

This step:

* Aligns expression matrices based on gene families
* Handles missing genes with NaN values
* Reports dimensions and recommendations

Step 6: Optimize Number of Core Components
-------------------------------------------

Use Cohen's d effect size metric to automatically determine the optimal number of core components.

.. code-block:: python

   # Optimize number of core components
   optimal_num_core_components = multiModulon.optimize_number_of_core_components(
       metric='effect_size',       # Use Cohen's d effect size
       step=5,                     # Test k = 5, 10, 15, 20, ...
       save_path=output_dir,       # Save plots
       fig_size=(7, 5),           # Figure size
   )

The optimization:

* Tests different numbers of core components
* Evaluates using Cohen's d effect size
* Selects optimal k based on interpretability
* Saves optimization plots

Step 7: Optimize Number of Unique Components
---------------------------------------------

Determine the optimal number of unique (species-specific) components for each strain.

.. code-block:: python

   # Optimize unique components for each species
   optimal_unique, optimal_total = multiModulon.optimize_number_of_unique_components(
       optimal_num_core_components=optimal_num_core_components,
       step=5,
       save_path=output_dir,
       fig_size=(7, 5)
   )

This process:

* Tests different numbers of unique components per species
* Evaluates component quality using effect size
* Returns optimal numbers for each species

Step 8: Run Robust Multi-view ICA
---------------------------------

Perform robust multi-view ICA with multiple runs and clustering to identify consistent components.

.. code-block:: python

   # Run robust multi-view ICA
   M_matrices, A_matrices = multiModulon.run_robust_multiview_ica(
       a=optimal_total,                 # Total components per species
       c=optimal_num_core_components,   # Number of core components
       num_runs=20,                     # Number of runs for robustness
       seed=42                          # Random seed
   )

The robust ICA:

* Runs ICA multiple times with different initializations
* Clusters components across runs
* Selects consistent components
* Generates final M (gene weights) and A (activities) matrices

Step 9: Optimize M Matrix Thresholds
------------------------------------

Calculate thresholds for binarizing the M matrices using Otsu's method.

.. code-block:: python

   # Optimize thresholds for each component
   multiModulon.optimize_M_thresholds(
       method="Otsu's method", 
       quantile_threshold=95
   )

This creates:

* Component-specific thresholds
* Binarized presence matrices
* Statistics on genes per component

Step 10: Save Results
---------------------

Save the complete MultiModulon object for future use.

.. code-block:: python

   # Save to compressed JSON format
   multiModulon.save_to_json_multimodulon("multiModulon_results.json.gz")
   
   # Load saved object
   multiModulon = MultiModulon.load_json_multimodulon("multiModulon_results.json.gz")

Summary
-------

This workflow covers the complete pipeline from raw data to interpretable multi-species regulatory modules:

1. Data loading and validation
2. Ortholog detection via BBH
3. Gene alignment across species
4. Expression matrix alignment
5. Component number optimization
6. Robust multi-view ICA
7. Threshold optimization
8. Results storage

The output includes:

* Core components conserved across species
* Unique components specific to each species
* Gene membership for each component
* Activity profiles across samples
* Optimized thresholds for interpretation