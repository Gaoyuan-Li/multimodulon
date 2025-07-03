Multi-view ICA
==============

This section covers running multi-view Independent Component Analysis (ICA) to identify core (conserved) and unique (species-specific) regulatory modules.

Overview
--------

Multi-view ICA decomposes expression matrices from multiple species into:

* **Core components** - Shared regulatory modules conserved across species
* **Unique components** - Species-specific regulatory modules

The analysis produces:

* **M matrices** - Gene weights (mixing matrices)
* **A matrices** - Sample activities (source signals)

Basic Multi-view ICA
--------------------

.. py:method:: MultiModulon.run_multiview_ica(**kwargs)

   Run multi-view ICA on aligned expression matrices.

   :param a: Number of components per species. Can be:
      - int: Same number for all species
      - Dict[str, int]: Different numbers per species
   :param int c: Number of core (shared) components
   :param str mode: 'gpu' or 'cpu' (default: 'gpu')
   :param float effect_size_threshold: Cohen's d threshold for filtering (optional)
   :param int num_top_gene: Number of top genes for Cohen's d calculation (default: 20)
   
   :return: Dictionary mapping species to M matrices
   :rtype: Dict[str, pd.DataFrame]

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   # Run with equal components per species
   M_matrices = mm.run_multiview_ica(
       a=50,      # 50 total components per species
       c=20,      # 20 core components
       mode='gpu'
   )
   
   # Run with different components per species
   M_matrices = mm.run_multiview_ica(
       a={'Species1': 50, 'Species2': 60, 'Species3': 45},
       c=20,
       mode='gpu'
   )

Component Filtering
~~~~~~~~~~~~~~~~~~~

Filter components by effect size:

.. code-block:: python

   # Only keep components with Cohen's d > 5
   M_matrices = mm.run_multiview_ica(
       a=60,
       c=25,
       effect_size_threshold=5,
       num_top_gene=20
   )
   
   # Check how many components passed the filter
   for species, M in M_matrices.items():
       print(f"{species}: {M.shape[1]} components retained")

Robust Multi-view ICA
---------------------

For more reliable results, use robust ICA with multiple runs:

.. py:method:: MultiModulon.run_robust_multiview_ica(**kwargs)

   Run multi-view ICA multiple times and cluster results for robustness.

   :param Dict[str, int] a: Total components per species
   :param int c: Number of core components
   :param int num_runs: Number of ICA runs (default: 100)
   :param str mode: 'gpu' or 'cpu' (default: 'gpu')
   :param int seed: Random seed (default: 42)
   :param float effect_size_threshold: Cohen's d threshold for all components
   :param float effect_size_threshold_core: Threshold for core components only
   :param float effect_size_threshold_unique: Threshold for unique components only
   :param int num_top_gene: Number of top genes for Cohen's d
   
   :return: Tuple of (M_matrices, A_matrices)
   :rtype: Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]

Robust ICA Example
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Run robust ICA with 100 runs
   M_matrices, A_matrices = mm.run_robust_multiview_ica(
       a={'Species1': 50, 'Species2': 60},
       c=20,
       num_runs=100,
       mode='gpu',
       seed=42
   )
   
   # Access results
   M_species1 = M_matrices['Species1']
   A_species1 = A_matrices['Species1']
   
   print(f"M matrix shape: {M_species1.shape}")
   print(f"A matrix shape: {A_species1.shape}")

Different Thresholds
~~~~~~~~~~~~~~~~~~~~

Apply different thresholds to core and unique components:

.. code-block:: python

   # Stricter threshold for core, looser for unique
   M_matrices, A_matrices = mm.run_robust_multiview_ica(
       a={'Species1': 50, 'Species2': 60},
       c=20,
       num_runs=100,
       effect_size_threshold_core=7,    # Strict for core
       effect_size_threshold_unique=3,  # Permissive for unique
       num_top_gene=20
   )

Understanding the Results
-------------------------

M Matrix (Gene Weights)
~~~~~~~~~~~~~~~~~~~~~~~

The M matrix contains gene weights for each component:

.. code-block:: python

   M = M_matrices['Species1']
   
   # Structure:
   # Rows: Genes (aligned across species)
   # Columns: Components (Core_1, Core_2, ..., Unique_1, ...)
   
   # Get top genes for a component
   component = 'Core_1'
   weights = M[component].sort_values(ascending=False)
   
   print(f"Top 10 genes in {component}:")
   print(weights.head(10))
   
   print(f"\nBottom 10 genes in {component}:")
   print(weights.tail(10))

A Matrix (Sample Activities)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The A matrix contains component activities across samples:

.. code-block:: python

   A = A_matrices['Species1']
   
   # Structure:
   # Rows: Components
   # Columns: Samples
   
   # Get activity profile for a component
   component = 'Core_1'
   activities = A.loc[component]
   
   # Find samples with high activity
   high_activity_samples = activities[activities > 5].index
   print(f"Samples with high {component} activity:")
   print(high_activity_samples.tolist())

Component Types
~~~~~~~~~~~~~~~

Components are labeled by type:

.. code-block:: python

   # List all components
   all_components = M.columns.tolist()
   
   # Separate by type
   core_components = [c for c in all_components if c.startswith('Core_')]
   unique_components = [c for c in all_components if c.startswith('Unique_')]
   
   print(f"Core components: {len(core_components)}")
   print(f"Unique components: {len(unique_components)}")

Generating Activity Matrices
----------------------------

After running ICA, generate A matrices from M and X:

.. py:method:: MultiModulon.generate_A()

   Generate A matrices (M.T @ X) for all species.
   
   **Example:**
   
   .. code-block:: python
      
      # Generate A matrices after ICA
      mm.generate_A()
      
      # Access generated matrices
      for species in mm.species:
          A = mm[species].A
          print(f"{species} activities: {A.shape}")

This is useful when:

* You've loaded pre-computed M matrices
* You want to recalculate activities after filtering

Advanced Usage
--------------

Custom Parameters
~~~~~~~~~~~~~~~~~

Fine-tune the ICA algorithm:

.. code-block:: python

   # Direct access to underlying function
   from multimodulon.multiview_ica import run_multiview_ica
   
   M_matrices = run_multiview_ica(
       species_X_matrices={s: mm[s].X for s in mm.species},
       a_values={'Species1': 50, 'Species2': 60},
       c=20,
       mode='gpu',
       max_iter=10000,      # More iterations
       learning_rate=0.01,  # Custom learning rate
       batch_size=None,     # Full batch
       seed=42
   )

GPU vs CPU Mode
~~~~~~~~~~~~~~~

Choose based on your system:

.. code-block:: python

   import torch
   
   # Check GPU availability
   if torch.cuda.is_available():
       print("GPU available - using GPU mode")
       mode = 'gpu'
   else:
       print("No GPU - using CPU mode")
       mode = 'cpu'
   
   # Run ICA
   M_matrices = mm.run_multiview_ica(
       a=50,
       c=20,
       mode=mode
   )

Handling Large Datasets
~~~~~~~~~~~~~~~~~~~~~~~

For memory efficiency:

.. code-block:: python

   # Use batched processing
   M_matrices = run_multiview_ica(
       species_X_matrices=X_matrices,
       a_values=a_values,
       c=20,
       batch_size=256,  # Process in batches
       mode='gpu'
   )

Quality Control
---------------

Assess ICA Results
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Calculate explained variance
   explained_var = mm.calculate_explained_variance()
   for species, var in explained_var.items():
       print(f"{species}: {var:.1%} variance explained")
   
   # Check component effect sizes
   from multimodulon.multiview_ica_optimization import calculate_average_effect_sizes
   
   effect_sizes = calculate_average_effect_sizes(
       M_matrices,
       num_top_gene=20
   )
   
   # Plot effect size distribution
   import matplotlib.pyplot as plt
   
   plt.hist(effect_sizes, bins=30)
   plt.xlabel("Cohen's d effect size")
   plt.ylabel("Number of components")
   plt.axvline(x=5, color='r', linestyle='--', label='Threshold')
   plt.legend()
   plt.show()

Component Correlation
~~~~~~~~~~~~~~~~~~~~~

Check independence of components:

.. code-block:: python

   # Within species
   M = M_matrices['Species1']
   corr_matrix = M.corr()
   
   # High correlation indicates redundancy
   import seaborn as sns
   
   plt.figure(figsize=(10, 8))
   sns.heatmap(corr_matrix, cmap='coolwarm', center=0)
   plt.title("Component correlation within Species1")
   plt.show()
   
   # Across species (for core components)
   core_comps = [c for c in M.columns if c.startswith('Core_')]
   M1_core = M_matrices['Species1'][core_comps]
   M2_core = M_matrices['Species2'][core_comps]
   
   # Compare matching components
   for comp in core_comps:
       corr = M1_core[comp].corr(M2_core[comp])
       print(f"{comp} correlation: {corr:.3f}")

Troubleshooting
---------------

**ICA doesn't converge:**

.. code-block:: python

   # Increase iterations
   M_matrices = mm.run_multiview_ica(
       a=50,
       c=20,
       max_iter=20000  # More iterations
   )
   
   # Or try different initialization
   M_matrices = mm.run_multiview_ica(
       a=50,
       c=20,
       seed=123  # Different random seed
   )

**Memory errors:**

.. code-block:: python

   # Reduce batch size or switch to CPU
   M_matrices = mm.run_multiview_ica(
       a=50,
       c=20,
       mode='cpu',
       batch_size=128  # Smaller batches
   )

**Poor component quality:**

.. code-block:: python

   # Try robust ICA with more runs
   M_matrices, A_matrices = mm.run_robust_multiview_ica(
       a={'Species1': 50, 'Species2': 60},
       c=20,
       num_runs=200,  # More runs
       effect_size_threshold=7  # Stricter threshold
   )

Best Practices
--------------

1. **Always use robust ICA** for final results (100+ runs)
2. **Check effect sizes** to ensure biological relevance
3. **Validate core components** across species
4. **Save intermediate results** for large analyses
5. **Use GPU mode** when available for speed

Next Steps
----------

After running ICA:

1. :doc:`visualization` - Visualize components
2. :doc:`examples/basic_workflow` - Complete workflow example
3. Biological interpretation - Analyze gene sets and activities