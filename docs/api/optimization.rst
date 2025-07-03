Optimization API Reference
==========================

This section provides detailed API documentation for optimization functions.

Core Component Optimization
---------------------------

.. autofunction:: multimodulon.optimization.optimize_number_of_core_components

   Detailed optimization of core component numbers using cross-validation.

   **Parameters:**
   
   * **multimodulon** (*MultiModulon*) -- MultiModulon instance with loaded data
   * **max_k** (*int*) -- Maximum number of core components to test
   * **step** (*int*) -- Step size for k candidates (default: 5)
   * **max_a_per_view** (*int*) -- Maximum components per view (default: 100)
   * **train_frac** (*float*) -- Fraction of data for training (default: 0.75)
   * **num_runs** (*int*) -- Number of cross-validation runs (default: 1)
   * **mode** (*str*) -- 'gpu' or 'cpu' (default: 'gpu')
   * **seed** (*int*) -- Random seed (default: 42)
   * **save_plot** (*str*) -- Path to save metric vs k plot (optional)
   * **metric** (*str*) -- 'nre' or 'effect_size' (default: 'nre')
   * **threshold** (*float*) -- Reconstruction threshold for NRE (optional)
   * **effect_size_threshold** (*float*) -- Cohen's d threshold (default: 5)
   * **num_top_gene** (*int*) -- Top genes for effect size (default: 20)
   * **save_path** (*str*) -- Directory to save plots (optional)
   * **fig_size** (*tuple*) -- Figure size (default: (5, 3))
   * **font_path** (*str*) -- Font file path (optional)
   
   **Returns:**
   
   * **optimal_k** (*int*) -- Optimal number of core components
   * **scores** (*Dict[int, float]*) -- Metric scores for each k tested
   
   **Example:**
   
   .. code-block:: python
      
      # NRE-based optimization
      optimal_k, scores = optimize_number_of_core_components(
          multimodulon=mm,
          max_k=40,
          step=5,
          metric='nre',
          num_runs=3,
          save_plot='nre_optimization.png'
      )
      
      # Effect size-based optimization
      optimal_k, scores = optimize_number_of_core_components(
          multimodulon=mm,
          max_k=40,
          step=5,
          metric='effect_size',
          effect_size_threshold=5,
          num_top_gene=20,
          save_plot='effect_size_optimization.png'
      )

Unique Component Optimization
-----------------------------

.. autofunction:: multimodulon.optimization.optimize_number_of_unique_components

   Optimize unique components for each species given fixed core components.

   **Parameters:**
   
   * **multimodulon** (*MultiModulon*) -- MultiModulon instance
   * **optimal_num_core_components** (*int*) -- Number of core components
   * **step** (*int*) -- Step size for testing (default: 5)
   * **mode** (*str*) -- 'gpu' or 'cpu' (default: 'gpu')
   * **seed** (*int*) -- Random seed (default: 42)
   * **save_plots** (*str*) -- Directory to save plots (optional)
   * **effect_size_threshold** (*float*) -- Cohen's d threshold (default: 5)
   * **num_top_gene** (*int*) -- Top genes for effect size (default: 20)
   * **save_path** (*str*) -- Directory path for plots (optional)
   * **fig_size** (*tuple*) -- Figure size (default: (5, 3))
   * **font_path** (*str*) -- Font file path (optional)
   
   **Returns:**
   
   * **optimal_unique** (*Dict[str, int]*) -- Optimal unique components per species
   * **optimal_total** (*Dict[str, int]*) -- Total components per species
   
   **Example:**
   
   .. code-block:: python
      
      optimal_unique, optimal_total = optimize_number_of_unique_components(
          multimodulon=mm,
          optimal_num_core_components=20,
          step=5,
          effect_size_threshold=5,
          save_plots='unique_optimization/'
      )
      
      # Results
      for species in mm.species:
          print(f"{species}: {optimal_unique[species]} unique, "
                f"{optimal_total[species]} total")

Multi-view ICA Optimization Module
----------------------------------

Effect Size Calculation
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: multimodulon.multiview_ica_optimization.calculate_cohens_d_effect_size

   Calculate Cohen's d effect size for a single component.

   **Parameters:**
   
   * **weight_vector** (*np.ndarray*) -- 1D array of gene weights
   * **seed** (*int*) -- Random seed (kept for compatibility)
   * **num_top_gene** (*int*) -- Number of top genes (default: 20)
   
   **Returns:**
   
   * **effect_size** (*float*) -- Cohen's d value
   
   **Example:**
   
   .. code-block:: python
      
      from multimodulon.multiview_ica_optimization import calculate_cohens_d_effect_size
      
      # For a single component
      weights = M_matrix.iloc[:, 0].values
      d = calculate_cohens_d_effect_size(weights, num_top_gene=20)
      print(f"Cohen's d = {d:.2f}")

.. autofunction:: multimodulon.multiview_ica_optimization.calculate_average_effect_sizes

   Calculate mean effect sizes across species for each component.

   **Parameters:**
   
   * **M_matrices** (*Dict[str, pd.DataFrame]*) -- M matrices per species
   * **seed** (*int*) -- Random seed (default: 42)
   * **num_top_gene** (*int*) -- Top genes (default: 20)
   
   **Returns:**
   
   * **effect_sizes** (*List[float]*) -- Mean effect size per component

NRE Calculation
~~~~~~~~~~~~~~~

.. autofunction:: multimodulon.multiview_ica_optimization.calculate_nre_proper

   Calculate Normalized Reconstruction Error for core components.

   **Parameters:**
   
   * **S_matrices** (*List[pd.DataFrame]*) -- Source signal matrices
   * **k_core** (*int*) -- Number of core components to evaluate
   
   **Returns:**
   
   * **nre** (*float*) -- NRE score (lower is better)

Optimization Runner
~~~~~~~~~~~~~~~~~~~

.. autofunction:: multimodulon.multiview_ica_optimization.run_nre_optimization

   Run optimization using specified metric with cross-validation.

   **Parameters:**
   
   * **species_X_matrices** (*Dict[str, pd.DataFrame]*) -- Expression matrices
   * **k_candidates** (*List[int]*) -- Core component numbers to test
   * **max_a_per_view** (*int*) -- Max components per species
   * **train_frac** (*float*) -- Training fraction (default: 0.75)
   * **num_runs** (*int*) -- CV runs (default: 1)
   * **mode** (*str*) -- 'gpu' or 'cpu' (default: 'gpu')
   * **seed** (*int*) -- Random seed (default: 42)
   * **metric** (*str*) -- 'nre' or 'effect_size' (default: 'nre')
   * **threshold** (*float*) -- NRE threshold (optional)
   * **effect_size_threshold** (*float*) -- Cohen's d threshold (default: 5)
   * **num_top_gene** (*int*) -- Top genes (default: 20)
   * **fig_size** (*tuple*) -- Figure size (default: (5, 3))
   * **font_path** (*str*) -- Font path (optional)
   
   **Returns:**
   
   * **optimal_k** (*int*) -- Optimal core components
   * **mean_scores** (*Dict[int, float]*) -- Mean scores per k
   * **all_scores** (*Dict[int, List]*) -- All scores per k
   * **std_scores** (*Dict[int, List]*) -- Standard deviations
   * **fig** (*plt.Figure*) -- Optimization plot

Usage Examples
--------------

Complete Optimization Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from multimodulon import MultiModulon
   from multimodulon.optimization import (
       optimize_number_of_core_components,
       optimize_number_of_unique_components
   )
   
   # Initialize and prepare data
   mm = MultiModulon("Input_Data")
   mm.generate_BBH()
   mm.align_genes()
   mm.generate_X("Output_Gene_Info")
   
   # Step 1: Optimize core components
   optimal_core, core_scores = optimize_number_of_core_components(
       multimodulon=mm,
       max_k=50,
       step=5,
       metric='effect_size',
       effect_size_threshold=5,
       num_runs=5,
       save_path='optimization_results/'
   )
   
   print(f"Optimal core components: {optimal_core}")
   
   # Step 2: Optimize unique components
   optimal_unique, optimal_total = optimize_number_of_unique_components(
       multimodulon=mm,
       optimal_num_core_components=optimal_core,
       step=5,
       effect_size_threshold=5,
       save_path='optimization_results/'
   )
   
   # Display results
   for species in mm.species:
       n_unique = optimal_unique[species]
       n_total = optimal_total[species]
       n_core = optimal_core
       print(f"{species}: {n_core} core + {n_unique} unique = {n_total} total")

Custom Metric Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def custom_metric(M_matrices, num_top_gene=20):
       """Custom metric for component evaluation."""
       scores = []
       
       for species, M in M_matrices.items():
           for col in M.columns:
               weights = M[col].values
               # Custom calculation
               score = your_custom_calculation(weights)
               scores.append(score)
       
       return np.mean(scores)
   
   # Use in optimization
   from multimodulon.multiview_ica_optimization import run_nre_optimization
   
   # Modify run_nre_optimization to use custom metric
   # or implement your own optimization loop

Analyzing Optimization Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   
   # Plot optimization landscape
   k_values = sorted(core_scores.keys())
   scores = [core_scores[k] for k in k_values]
   
   plt.figure(figsize=(8, 5))
   plt.plot(k_values, scores, 'o-', linewidth=2, markersize=8)
   plt.xlabel('Number of core components')
   plt.ylabel('Metric score')
   plt.title('Core component optimization')
   plt.grid(True, alpha=0.3)
   
   # Mark optimal
   optimal_score = core_scores[optimal_core]
   plt.axvline(optimal_core, color='red', linestyle='--', label=f'Optimal k={optimal_core}')
   plt.legend()
   plt.tight_layout()
   plt.savefig('optimization_landscape.png', dpi=300)
   
   # Analyze stability
   if isinstance(scores[0], list):  # Multiple runs
       means = [np.mean(s) for s in scores]
       stds = [np.std(s) for s in scores]
       
       plt.figure(figsize=(8, 5))
       plt.errorbar(k_values, means, yerr=stds, fmt='o-', capsize=5)
       plt.xlabel('Number of core components')
       plt.ylabel('Metric score (mean ± std)')
       plt.title('Optimization stability across runs')
       plt.tight_layout()
       plt.savefig('optimization_stability.png', dpi=300)