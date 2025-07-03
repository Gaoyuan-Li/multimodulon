Optimization of Dimensions
==========================

This section covers the optimization of component numbers for multi-view ICA, including both core (shared) and unique (species-specific) components.

Overview
--------

Choosing the right number of components is crucial for meaningful results. MultiModulon provides automated optimization methods to determine:

1. **Optimal number of core components** - Shared across all species
2. **Optimal number of unique components** - Specific to each species

Two optimization metrics are available:

* **NRE (Normalized Reconstruction Error)** - Measures reconstruction quality
* **Cohen's d Effect Size** - Measures biological interpretability

Optimizing Core Components
--------------------------

.. py:method:: MultiModulon.optimize_number_of_core_components(**kwargs)

   Optimize the number of core (shared) components across species.

   :param int max_k: Maximum number of core components to test
   :param int step: Step size for k candidates (default: 5)
   :param int max_a_per_view: Maximum components per species (default: 100)
   :param float train_frac: Fraction of data for training (default: 0.75)
   :param int num_runs: Number of cross-validation runs (default: 1)
   :param str mode: Computation mode 'gpu' or 'cpu' (default: 'gpu')
   :param int seed: Random seed for reproducibility (default: 42)
   :param str metric: Optimization metric 'nre' or 'effect_size' (default: 'nre')
   :param float effect_size_threshold: Cohen's d threshold (default: 5)
   :param int num_top_gene: Number of top genes for Cohen's d (default: 20)
   :param str save_path: Directory to save optimization plot
   :param tuple fig_size: Figure size as (width, height) (default: (5, 3))
   :param str font_path: Path to font file for plots
   
   :return: Tuple of (optimal_num_core_components, metric_scores)
   :rtype: Tuple[int, Dict[int, float]]

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   # Optimize using NRE metric
   optimal_core, scores = mm.optimize_number_of_core_components(
       max_k=30,
       step=5,
       mode='gpu',
       save_plot="core_optimization.png"
   )
   
   print(f"Optimal number of core components: {optimal_core}")

Using Effect Size Metric
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Optimize using Cohen's d effect size
   optimal_core, scores = mm.optimize_number_of_core_components(
       max_k=30,
       step=5,
       metric='effect_size',
       effect_size_threshold=5,  # Minimum Cohen's d
       num_top_gene=20,          # Top genes to consider
       save_plot="effect_size_optimization.png"
   )

Understanding the Metrics
~~~~~~~~~~~~~~~~~~~~~~~~~

**NRE (Normalized Reconstruction Error):**

* Measures how well core components reconstruct the data
* Lower values are better
* Good for capturing maximum variance
* May include noise components

**Cohen's d Effect Size:**

* Measures separation between top genes and others
* Higher values indicate more specific components
* Better for biological interpretability
* Filters out non-specific components

Cross-validation
~~~~~~~~~~~~~~~~

For robust estimates, use multiple runs:

.. code-block:: python

   # Multiple cross-validation runs
   optimal_core, scores = mm.optimize_number_of_core_components(
       max_k=30,
       step=5,
       train_frac=0.75,  # 75% train, 25% test
       num_runs=5,       # 5-fold cross-validation
       save_plot="cv_optimization.png"
   )
   
   # Scores now contain mean and std across runs
   for k, score in scores.items():
       print(f"k={k}: {score:.3f}")

Optimizing Unique Components
----------------------------

After determining core components, optimize unique components:

.. py:method:: MultiModulon.optimize_number_of_unique_components(**kwargs)

   Optimize the number of unique components for each species.

   :param int optimal_num_core_components: Number of core components (from previous step)
   :param int step: Step size for testing unique components (default: 5)
   :param str mode: Computation mode 'gpu' or 'cpu' (default: 'gpu')
   :param int seed: Random seed (default: 42)
   :param float effect_size_threshold: Cohen's d threshold (default: 5)
   :param int num_top_gene: Number of top genes for Cohen's d (default: 20)
   :param str save_path: Directory to save plots for each species
   :param tuple fig_size: Figure size (default: (5, 3))
   :param str font_path: Path to font file
   
   :return: Tuple of (optimal_unique_components, optimal_total_components)
   :rtype: Tuple[Dict[str, int], Dict[str, int]]

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   # Optimize unique components
   optimal_unique, optimal_total = mm.optimize_number_of_unique_components(
       optimal_num_core_components=20,  # From previous step
       step=5,
       save_plots="unique_optimization/"
   )
   
   # Results
   print("Optimal unique components per species:")
   for species, n_unique in optimal_unique.items():
       n_total = optimal_total[species]
       print(f"{species}: {n_unique} unique, {n_total} total")

How It Works
~~~~~~~~~~~~

For each species:

1. Tests different numbers of unique components
2. Runs ICA with fixed core + varying unique
3. Calculates mean Cohen's d for unique components
4. Selects number that maximizes interpretable components

Custom Thresholds
~~~~~~~~~~~~~~~~~

Different species may need different thresholds:

.. code-block:: python

   # Strict threshold for well-studied species
   optimal_unique_strict, _ = mm.optimize_number_of_unique_components(
       optimal_num_core_components=20,
       effect_size_threshold=7,  # Higher threshold
       save_plots="strict_optimization/"
   )
   
   # Permissive threshold for novel species  
   optimal_unique_permissive, _ = mm.optimize_number_of_unique_components(
       optimal_num_core_components=20,
       effect_size_threshold=3,  # Lower threshold
       save_plots="permissive_optimization/"
   )

Complete Optimization Workflow
------------------------------

Here's a complete optimization workflow:

.. code-block:: python

   # Step 1: Optimize core components
   print("Optimizing core components...")
   optimal_core, core_scores = mm.optimize_number_of_core_components(
       max_k=40,
       step=5,
       metric='effect_size',
       effect_size_threshold=5,
       num_runs=3,
       save_path="optimization_results/",
       fig_size=(6, 4)
   )
   print(f"Optimal core components: {optimal_core}")
   
   # Step 2: Optimize unique components
   print("\nOptimizing unique components...")
   optimal_unique, optimal_total = mm.optimize_number_of_unique_components(
       optimal_num_core_components=optimal_core,
       step=5,
       effect_size_threshold=5,
       save_path="optimization_results/",
       fig_size=(6, 4)
   )
   
   print("\nOptimization complete!")
   print(f"Core components: {optimal_core}")
   for species in mm.species:
       print(f"{species}: {optimal_unique[species]} unique, "
             f"{optimal_total[species]} total")
   
   # Step 3: Run ICA with optimal parameters
   print("\nRunning multi-view ICA with optimal parameters...")
   M_matrices, A_matrices = mm.run_robust_multiview_ica(
       a=optimal_total,
       c=optimal_core,
       num_runs=100,
       mode='gpu'
   )

Advanced Options
----------------

Manual Component Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~

Override automatic selection:

.. code-block:: python

   # Examine optimization curves
   import matplotlib.pyplot as plt
   
   # Plot scores
   k_values = list(scores.keys())
   score_values = list(scores.values())
   
   plt.plot(k_values, score_values, 'o-')
   plt.xlabel('Number of core components')
   plt.ylabel('Metric score')
   plt.show()
   
   # Manually select based on curve
   manual_core = 25  # Your choice

Effect Size Calculation
~~~~~~~~~~~~~~~~~~~~~~~

Understand how Cohen's d is calculated:

.. code-block:: python

   from multimodulon.multiview_ica_optimization import calculate_cohens_d_effect_size
   
   # For a single component
   gene_weights = M_matrix.iloc[:, 0].values  # First component
   effect_size = calculate_cohens_d_effect_size(
       gene_weights,
       num_top_gene=20
   )
   print(f"Cohen's d = {effect_size:.2f}")

Parallel Optimization
~~~~~~~~~~~~~~~~~~~~~

Speed up optimization with parallel runs:

.. code-block:: python

   # GPU mode automatically parallelizes
   # For CPU mode, consider multiprocessing:
   
   from multiprocessing import Pool
   
   def optimize_k(k):
       # Run optimization for single k
       return k, run_ica_and_calculate_metric(k)
   
   with Pool(processes=4) as pool:
       results = pool.map(optimize_k, range(5, 35, 5))

Troubleshooting
---------------

**Optimization takes too long:**

.. code-block:: python

   # Reduce search space
   optimal_core, _ = mm.optimize_number_of_core_components(
       max_k=20,      # Lower maximum
       step=10,       # Larger steps
       num_runs=1,    # Single run
       mode='gpu'     # Use GPU
   )

**No clear optimum:**

.. code-block:: python

   # Try different metric
   # If NRE keeps decreasing, try effect_size
   # If effect_size is flat, try different threshold
   
   # Save all plots to compare
   for metric in ['nre', 'effect_size']:
       for threshold in [3, 5, 7]:
           optimal, scores = mm.optimize_number_of_core_components(
               metric=metric,
               effect_size_threshold=threshold,
               save_plot=f"opt_{metric}_t{threshold}.png"
           )

**Different results between runs:**

.. code-block:: python

   # Increase number of runs for stability
   optimal_core, scores = mm.optimize_number_of_core_components(
       num_runs=10,     # More runs
       seed=42,         # Fix seed
       train_frac=0.8   # More training data
   )

Best Practices
--------------

1. **Start with effect_size metric** - More biologically relevant
2. **Use multiple runs** - At least 3-5 for reliability  
3. **Inspect plots** - Don't just trust automatic selection
4. **Consider biology** - Known pathways can guide selection
5. **Validate results** - Check if components make biological sense

Next Steps
----------

After optimization:

1. :doc:`multiview_ica` - Run ICA with optimal parameters
2. :doc:`visualization` - Visualize and interpret components
3. :doc:`examples/optimization_tutorial` - Detailed optimization tutorial