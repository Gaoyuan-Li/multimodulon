Optimization Tutorial
=====================

This tutorial covers the optimization strategies for determining the optimal number of core and unique components.

Understanding Component Optimization
------------------------------------

The quality of multi-view ICA results depends critically on choosing the right number of components. Too few components merge distinct regulatory modules, while too many create noise and redundancy.

Optimization Metrics
--------------------

Normalized Reconstruction Error (NRE)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NRE measures how well the components reconstruct the original data:

.. code-block:: python

   # Optimize using NRE metric
   optimal_core, scores = multiModulon.optimize_number_of_core_components(
       metric='nre',
       max_k=50,
       step=5,
       train_frac=0.75,  # 75% train, 25% test
       num_runs=3,       # Cross-validation runs
       mode='gpu'
   )

Advantages:
* Captures maximum variance
* Well-established metric
* Good for initial exploration

Limitations:
* May include noise components
* Doesn't consider biological interpretability

Cohen's d Effect Size
~~~~~~~~~~~~~~~~~~~~~

Cohen's d measures the separation between top genes and the rest:

.. code-block:: python

   # Optimize using effect size metric
   optimal_core, scores = multiModulon.optimize_number_of_core_components(
       metric='effect_size',
       max_k=50,
       step=5,
       effect_size_threshold=5,    # Minimum Cohen's d
       num_top_gene=20,           # Top genes to consider
       num_runs=3,
       save_path='optimization_plots/'
   )

Advantages:
* Focuses on interpretable components
* Filters noise automatically
* Better biological relevance

Core Component Optimization
---------------------------

Step-by-Step Process
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Step 1: Coarse search
   coarse_optimal, coarse_scores = multiModulon.optimize_number_of_core_components(
       max_k=60,
       step=10,  # Large steps: 10, 20, 30, ...
       metric='effect_size',
       effect_size_threshold=5,
       save_path='coarse_search.png'
   )
   
   # Step 2: Fine search around optimum
   fine_optimal, fine_scores = multiModulon.optimize_number_of_core_components(
       max_k=coarse_optimal + 10,
       min_k=max(5, coarse_optimal - 10),
       step=2,  # Small steps
       metric='effect_size',
       effect_size_threshold=5,
       save_path='fine_search.png'
   )

Interpreting Results
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Analyze optimization curve
   import matplotlib.pyplot as plt
   
   k_values = sorted(scores.keys())
   score_values = [scores[k] for k in k_values]
   
   plt.figure(figsize=(8, 5))
   plt.plot(k_values, score_values, 'o-', linewidth=2)
   plt.xlabel('Number of core components')
   plt.ylabel('Average components above threshold')
   plt.title('Core Component Optimization')
   plt.grid(True, alpha=0.3)
   
   # Mark optimal
   plt.axvline(optimal_core, color='red', linestyle='--', 
               label=f'Optimal k={optimal_core}')
   plt.legend()
   plt.show()

Look for:
* Clear peak or plateau
* Stability across nearby values
* Balance between too few and too many

Unique Component Optimization
-----------------------------

Basic Optimization
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Optimize unique components
   optimal_unique, optimal_total = multiModulon.optimize_number_of_unique_components(
       optimal_num_core_components=optimal_core,
       step=5,
       effect_size_threshold=1,  # Can be lower than core
       save_path='unique_optimization/'
   )
   
   # Review results
   for species in multiModulon.species:
       print(f"{species}:")
       print(f"  Core: {optimal_core}")
       print(f"  Unique: {optimal_unique[species]}")
       print(f"  Total: {optimal_total[species]}")

Species-Specific Thresholds
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Different species may need different criteria:

.. code-block:: python

   # Custom thresholds per species type
   threshold_map = {
       'well_studied': 5,      # Higher for well-annotated species
       'novel': 3,             # Lower for less-studied species
       'environmental': 4      # Medium for environmental isolates
   }
   
   # Run optimization with custom thresholds
   species_results = {}
   for species, category in species_categories.items():
       result = optimize_unique_for_species(
           species=species,
           core_k=optimal_core,
           threshold=threshold_map[category]
       )
       species_results[species] = result

Advanced Optimization Strategies
--------------------------------

Multi-Metric Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~

Combine multiple metrics for robust selection:

.. code-block:: python

   # Run optimization with both metrics
   nre_optimal, nre_scores = multiModulon.optimize_number_of_core_components(
       metric='nre',
       max_k=50,
       step=5
   )
   
   effect_optimal, effect_scores = multiModulon.optimize_number_of_core_components(
       metric='effect_size',
       max_k=50,
       step=5,
       effect_size_threshold=5
   )
   
   # Find consensus
   if abs(nre_optimal - effect_optimal) <= 5:
       # Close agreement - use effect size result
       final_optimal = effect_optimal
   else:
       # Disagreement - need manual inspection
       print(f"NRE suggests k={nre_optimal}")
       print(f"Effect size suggests k={effect_optimal}")
       # Examine both results before deciding

Stability Analysis
~~~~~~~~~~~~~~~~~~

Test optimization stability:

.. code-block:: python

   # Multiple optimization runs
   n_iterations = 5
   optimal_values = []
   
   for i in range(n_iterations):
       optimal, _ = multiModulon.optimize_number_of_core_components(
           metric='effect_size',
           max_k=50,
           step=5,
           seed=i * 42,  # Different seeds
           num_runs=3
       )
       optimal_values.append(optimal)
   
   # Check consistency
   import numpy as np
   print(f"Optimal values: {optimal_values}")
   print(f"Mean: {np.mean(optimal_values):.1f}")
   print(f"Std: {np.std(optimal_values):.1f}")
   
   # If std is low, optimization is stable

Biological Validation
~~~~~~~~~~~~~~~~~~~~~

Validate optimization using known biology:

.. code-block:: python

   # Test different k values around optimum
   test_k_values = [optimal_core - 5, optimal_core, optimal_core + 5]
   
   for k in test_k_values:
       print(f"\nTesting k={k}:")
       
       # Run ICA
       M_matrices = multiModulon.run_multiview_ica(
           a=60,  # Fixed total
           c=k,   # Variable core
           effect_size_threshold=5
       )
       
       # Check known regulons
       known_regulons = ['fur', 'crp', 'fnr', 'arcA']
       for regulon in known_regulons:
           # Search for regulon in components
           found = search_regulon_in_components(M_matrices, regulon)
           print(f"  {regulon}: {'Found' if found else 'Not found'}")

Troubleshooting Optimization
----------------------------

No Clear Optimum
~~~~~~~~~~~~~~~~~

If optimization curve is flat or monotonic:

.. code-block:: python

   # Try different parameters
   param_combinations = [
       {'effect_size_threshold': 3, 'num_top_gene': 20},
       {'effect_size_threshold': 5, 'num_top_gene': 20},
       {'effect_size_threshold': 5, 'num_top_gene': 30},
       {'effect_size_threshold': 7, 'num_top_gene': 20},
   ]
   
   results = {}
   for i, params in enumerate(param_combinations):
       optimal, scores = multiModulon.optimize_number_of_core_components(
           metric='effect_size',
           max_k=50,
           step=5,
           **params
       )
       results[f"params_{i}"] = (optimal, scores, params)
   
   # Compare results
   for key, (opt, _, params) in results.items():
       print(f"{key}: k={opt}, params={params}")

Computational Constraints
~~~~~~~~~~~~~~~~~~~~~~~~~

For large datasets:

.. code-block:: python

   # Reduce computational load
   # Option 1: Larger steps
   optimal_quick = multiModulon.optimize_number_of_core_components(
       max_k=60,
       step=15,  # Test fewer values
       num_runs=1,  # Single run
       mode='gpu'
   )
   
   # Option 2: Subsample data
   # Create smaller dataset for optimization
   subsample_fraction = 0.5
   # Run optimization on subset
   # Then validate on full data

Best Practices
--------------

1. **Start with effect_size metric** - More biologically relevant
2. **Use multiple runs** - Ensures stability (num_runs ≥ 3)
3. **Validate with biology** - Check if known modules are captured
4. **Document choices** - Record parameters and rationale
5. **Consider species differences** - May need different thresholds

Optimization Workflow Summary
-----------------------------

.. code-block:: python

   # Complete optimization workflow
   
   # 1. Optimize core components
   optimal_core, core_scores = multiModulon.optimize_number_of_core_components(
       metric='effect_size',
       max_k=60,
       step=5,
       effect_size_threshold=5,
       num_runs=5,
       save_path='optimization/core/'
   )
   
   print(f"Optimal core components: {optimal_core}")
   
   # 2. Optimize unique components
   optimal_unique, optimal_total = multiModulon.optimize_number_of_unique_components(
       optimal_num_core_components=optimal_core,
       step=5,
       effect_size_threshold=3,  # Lower threshold for unique
       save_path='optimization/unique/'
   )
   
   # 3. Display final configuration
   print("\nFinal component configuration:")
   print(f"Core components: {optimal_core}")
   for species in multiModulon.species:
       print(f"{species}: {optimal_unique[species]} unique, "
             f"{optimal_total[species]} total")
   
   # 4. Run ICA with optimal parameters
   M_matrices, A_matrices = multiModulon.run_robust_multiview_ica(
       a=optimal_total,
       c=optimal_core,
       num_runs=100,  # Many runs for final result
       seed=42
   )
   
   print("\nOptimization complete!")