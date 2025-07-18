Optimization of Dimensions
==========================

This section covers the optimization of component numbers for multi-view ICA, including both core (shared) and unique (species-specific) components.

Overview
--------

Choosing the right number of components is crucial for meaningful results. MultiModulon provides automated optimization methods to determine:

1. **Optimal number of core components** - Shared across all species
2. **Optimal number of unique components** - Specific to each species

Two optimization metrics are available:

* **Cohen's d Effect Size** - between top genes and others (default, and recommended)
* **NRE (Normalized Reconstruction Error)** - From paper https://proceedings.mlr.press/v216/pandeva23a.html)

Optimizing Core Components
--------------------------

.. py:method:: MultiModulon.optimize_number_of_core_components(**kwargs)

   Optimize the number of core (shared) components across species.

   :param int max_k: Maximum number of core components to test (Auto-determined)
   :param int step: Step size for k candidates (default: 5)
   :param int max_a_per_view: Maximum components per species (default: max_k)
   :param float train_frac: Fraction of data for training (default: 0.75)
   :param int num_runs: Number of cross-validation runs (default: 1)
   :param str mode: Computation mode 'gpu' or 'cpu' (default: 'gpu')
   :param int seed: Random seed for reproducibility (default: 42)
   :param str metric: Optimization metric 'nre' or 'effect_size' (default: 'effect_size')
   :param str save_path: Directory to save optimization plot
   :param tuple fig_size: Figure size as (width, height) (default: (5, 3))
   :param str font_path: Path to font file for plots
   
   :return: Tuple of (optimal_num_core_components, metric_scores)
   :rtype: Tuple[int, Dict[int, float]]

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   # Optimize using Cohen's d effect size
   optimal_core, scores = multiModulon.optimize_number_of_core_components(
       max_k=30,
       step=5,
       metric='effect_size',
       save_plot="effect_size_optimization.png"
   )

   print(f"Optimal number of core components: {optimal_core}")

Understanding the Metrics
~~~~~~~~~~~~~~~~~~~~~~~~~

**NRE (Normalized Reconstruction Error):**

* Measures how well core components reconstruct the data
* Lower values are better
* May include noise components

**Cohen's d Effect Size:**

* Measures separation between top genes and others
* Higher values indicate components with a more clear gene membership
* Better for biological interpretability
* Filters out noise components

Optimizing Unique Components
----------------------------

After determining core components, optimize unique components:

.. py:method:: MultiModulon.optimize_number_of_unique_components(**kwargs)

   Optimize the number of unique components for each species.

   :param int optimal_num_core_components: Number of core components (from previous step)
   :param int step: Step size for testing unique components (default: 5)
   :param str mode: Computation mode 'gpu' or 'cpu' (default: 'gpu')
   :param int seed: Random seed (default: 42)
   :param str save_path: Directory to save plots for each species
   :param tuple fig_size: Figure size (default: (5, 3))
   :param str font_path: Path to font file
   
   :return: Tuple of (optimal_unique_components, optimal_total_components)
   :rtype: Tuple[Dict[str, int], Dict[str, int]]

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   # Optimize unique components
   optimal_unique, optimal_total = multiModulon.optimize_number_of_unique_components(
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


Complete Optimization Workflow
------------------------------

Here's a complete optimization workflow:

.. code-block:: python

   # Step 1: Optimize core components
   print("Optimizing core components...")
   optimal_core, core_scores = multiModulon.optimize_number_of_core_components(
       max_k=40,
       step=5,
       metric='effect_size',
       num_runs=3,
       save_path="optimization_results/",
       fig_size=(6, 4)
   )
   print(f"Optimal core components: {optimal_core}")
   
   # Step 2: Optimize unique components
   print("\nOptimizing unique components...")
   optimal_unique, optimal_total = multiModulon.optimize_number_of_unique_components(
       optimal_num_core_components=optimal_core,
       step=5,
       save_path="optimization_results/",
       fig_size=(6, 4)
   )
   
   print("\nOptimization complete!")
   print(f"Core components: {optimal_core}")
   for species in multiModulon.species:
       print(f"{species}: {optimal_unique[species]} unique, "
             f"{optimal_total[species]} total")
   
   # Step 3: Run ICA with optimal parameters
   print("\nRunning multi-view ICA with optimal parameters...")
   M_matrices, A_matrices = multiModulon.run_robust_multiview_ica(
       a=optimal_total,
       c=optimal_core,
       num_runs=100,
       mode='gpu'
   )

Best Practices
--------------

1. **Start with effect_size metric** - More biologically relevant
2. **Use multiple runs** - At least 3-5 for reliability  
3. **Inspect plots** - Don't just trust automatic selection
4. **Validate results** - Check if components make biological sense

Next Steps
----------

After optimization:

1. :doc:`multiview_ica` - Run ICA with optimal parameters
2. :doc:`visualization` - Visualize and interpret components