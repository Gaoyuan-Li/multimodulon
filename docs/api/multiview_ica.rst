Multi-view ICA API Reference
============================

This section provides detailed API documentation for multi-view ICA functions.

Main ICA Function
-----------------

.. autofunction:: multimodulon.multiview_ica.run_multiview_ica

   Run multi-view ICA using PyTorch implementation.

   **Parameters:**
   
   * **species_X_matrices** (*Dict[str, pd.DataFrame]*) -- Aligned expression matrices
   * **a_values** (*Dict[str, int]*) -- Total components per species
   * **c** (*int*) -- Number of core (shared) components
   * **mode** (*str*) -- 'gpu' or 'cpu' (default: 'gpu')
   * **return_unmixing_matrices** (*bool*) -- Return W matrices (default: False)
   * **effect_size_threshold** (*float*) -- Filter components by Cohen's d (optional)
   * **num_top_gene** (*int*) -- Top genes for effect size (default: 20)
   * **max_iter** (*int*) -- Maximum iterations (default: 10000)
   * **learning_rate** (*float*) -- Learning rate (default: 0.01)
   * **batch_size** (*int*) -- Batch size (default: None - full batch)
   * **seed** (*int*) -- Random seed (default: 0)
   
   **Returns:**
   
   * **M_matrices** (*Dict[str, pd.DataFrame]*) -- Mixing matrices (gene weights)
   * **W_matrices** (*Dict[str, pd.DataFrame]*) -- Unmixing matrices (if requested)
   
   **Example:**
   
   .. code-block:: python
      
      from multimodulon.multiview_ica import run_multiview_ica
      
      # Prepare expression matrices
      X_matrices = {
          'Species1': species1_X,
          'Species2': species2_X
      }
      
      # Run ICA
      M_matrices = run_multiview_ica(
          species_X_matrices=X_matrices,
          a_values={'Species1': 50, 'Species2': 60},
          c=20,
          mode='gpu',
          effect_size_threshold=5
      )

Core ICA Implementation
-----------------------

.. autofunction:: multimodulon.multiview_ica.run_multi_view_ICA_on_datasets

   Low-level multi-view ICA implementation.

   **Parameters:**
   
   * **datasets** (*List[pd.DataFrame]*) -- List of expression matrices
   * **a_values** (*List[int]*) -- Components per view
   * **c** (*int*) -- Core components
   * **batch_size** (*int*) -- Batch size (optional)
   * **max_iter** (*int*) -- Max iterations (default: 10000)
   * **seed** (*int*) -- Random seed (default: 0)
   * **mode** (*str*) -- 'gpu' or 'cpu' (default: 'gpu')
   * **return_unmixing_matrices** (*bool*) -- Return W matrices (default: False)
   
   **Returns:**
   
   * **results** (*List[pd.DataFrame]*) -- Source signals or (S, W) tuples
   
   **Algorithm:**
   
   1. Whitens input data
   2. Optimizes orthogonal unmixing matrices
   3. Extracts independent components
   4. Orders components across views

MSIICA Model
------------

.. autoclass:: multimodulon.multiview_ica.MSIICA
   :members:
   :undoc-members:
   :show-inheritance:
   
   Multi-Subject Independent ICA model with orthogonal constraints.
   
   **Constructor Parameters:**
   
   * **n_in** (*List[int]*) -- Input dimensions per view
   * **n_out** (*List[int]*) -- Output dimensions per view
   * **U** (*torch.Tensor*) -- Initial unmixing matrices (optional)
   * **ortho** (*bool*) -- Use orthogonal constraint (default: True)
   
   **Methods:**
   
   .. automethod:: forward
   .. automethod:: loss
   
   **Example:**
   
   .. code-block:: python
      
      import torch
      from multimodulon.multiview_ica import MSIICA
      
      # Create model
      model = MSIICA(
          n_in=[100, 100],    # 100 components after whitening
          n_out=[50, 60],     # 50 and 60 total components
          ortho=True
      )
      
      # Forward pass
      Xw = [torch.randn(1, 100, 1000), torch.randn(1, 100, 1000)]
      S = model(Xw)

Helper Functions
----------------

Whitening
~~~~~~~~~

.. autofunction:: multimodulon.multiview_ica.whiten

   GPU-accelerated data whitening.

   **Parameters:**
   
   * **X** (*torch.Tensor*) -- Input data (1, features, samples)
   * **rank** (*int*) -- Number of components to retain
   
   **Returns:**
   
   * **K** (*torch.Tensor*) -- Whitening matrix
   * **Xw** (*torch.Tensor*) -- Whitened data
   
   **Example:**
   
   .. code-block:: python
      
      import torch
      from multimodulon.multiview_ica import whiten
      
      # Prepare data
      X = torch.randn(1, 3000, 500)  # 3000 genes, 500 samples
      
      # Whiten to 100 components
      K, Xw = whiten(X, rank=100)
      print(f"Whitened shape: {Xw.shape}")  # (1, 100, 500)

Component Ordering
~~~~~~~~~~~~~~~~~~

.. autofunction:: multimodulon.multiview_ica.find_ordering

   Find optimal ordering of components across views.

   **Parameters:**
   
   * **S_list** (*List[np.ndarray]*) -- Source signals per view
   
   **Returns:**
   
   * **u** (*int*) -- Number of unique components per view
   * **orders** (*List[np.ndarray]*) -- Component ordering
   * **vals** (*np.ndarray*) -- Correlation values
   
   **Algorithm:**
   
   Uses Hungarian algorithm to match components based on correlation.

Usage Examples
--------------

Basic Multi-view ICA
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from multimodulon import MultiModulon
   from multimodulon.multiview_ica import run_multiview_ica
   
   # Load and prepare data
   mm = MultiModulon("Input_Data")
   mm.generate_X("Output_Gene_Info")
   
   # Get aligned expression matrices
   X_matrices = {}
   for species in mm.species:
       X_matrices[species] = mm[species].X
   
   # Run multi-view ICA
   M_matrices = run_multiview_ica(
       species_X_matrices=X_matrices,
       a_values={'Species1': 50, 'Species2': 60, 'Species3': 55},
       c=25,  # 25 core components
       mode='gpu',
       max_iter=10000,
       seed=42
   )
   
   # Save results
   for species, M in M_matrices.items():
       M.to_csv(f"{species}_M_matrix.csv")

Component Filtering
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Filter by effect size
   M_filtered = run_multiview_ica(
       species_X_matrices=X_matrices,
       a_values={'Species1': 60, 'Species2': 70},
       c=30,
       effect_size_threshold=5,  # Keep only components with d > 5
       num_top_gene=20
   )
   
   # Check filtering results
   for species, M in M_filtered.items():
       n_original = 60 if species == 'Species1' else 70
       n_filtered = M.shape[1]
       print(f"{species}: {n_filtered}/{n_original} components retained")

GPU Memory Management
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import torch
   
   # Check GPU memory
   if torch.cuda.is_available():
       print(f"GPU: {torch.cuda.get_device_name()}")
       print(f"Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
   
   # Run with batching for large datasets
   M_matrices = run_multiview_ica(
       species_X_matrices=X_matrices,
       a_values={'Species1': 100, 'Species2': 120},
       c=50,
       mode='gpu',
       batch_size=256  # Process in batches
   )
   
   # Or use CPU for very large datasets
   M_matrices_cpu = run_multiview_ica(
       species_X_matrices=X_matrices,
       a_values={'Species1': 100, 'Species2': 120},
       c=50,
       mode='cpu'
   )

Custom ICA Parameters
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Fine-tune optimization
   M_matrices = run_multiview_ica(
       species_X_matrices=X_matrices,
       a_values={'Species1': 50, 'Species2': 60},
       c=20,
       mode='gpu',
       max_iter=20000,      # More iterations
       learning_rate=0.005,  # Slower learning
       seed=123
   )
   
   # Get unmixing matrices too
   M_matrices, W_matrices = run_multiview_ica(
       species_X_matrices=X_matrices,
       a_values={'Species1': 50, 'Species2': 60},
       c=20,
       return_unmixing_matrices=True
   )
   
   # Verify relationship: S = W @ X
   X1 = X_matrices['Species1']
   W1 = W_matrices['Species1']
   S1_computed = W1 @ X1
   print(f"Source signal shape: {S1_computed.shape}")

Robust ICA with Multiple Runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from multimodulon.multiview_ica import run_multiview_ica
   import numpy as np
   
   # Run ICA multiple times
   n_runs = 10
   all_M_matrices = []
   
   for run in range(n_runs):
       print(f"Run {run + 1}/{n_runs}")
       M_matrices = run_multiview_ica(
           species_X_matrices=X_matrices,
           a_values={'Species1': 50, 'Species2': 60},
           c=20,
           seed=run  # Different seed each run
       )
       all_M_matrices.append(M_matrices)
   
   # Analyze consistency
   # For each species, cluster components across runs
   from sklearn.cluster import AgglomerativeClustering
   
   for species in X_matrices.keys():
       # Stack all M matrices
       all_M = np.hstack([
           run[species].values for run in all_M_matrices
       ])
       
       # Cluster components
       clustering = AgglomerativeClustering(
           n_clusters=all_M_matrices[0][species].shape[1],
           linkage='average'
       )
       labels = clustering.fit_predict(all_M.T)
       
       print(f"{species}: Component clustering consistency")

Troubleshooting
---------------

**CUDA out of memory:**

.. code-block:: python

   # Clear GPU cache
   torch.cuda.empty_cache()
   
   # Use smaller batch size
   M = run_multiview_ica(
       species_X_matrices=X_matrices,
       a_values=a_values,
       c=c,
       batch_size=128  # Smaller batches
   )

**Convergence issues:**

.. code-block:: python

   # Monitor loss during training
   import matplotlib.pyplot as plt
   
   # Modify run_multi_view_ICA_on_datasets to return loss history
   # Then plot:
   plt.plot(loss_history)
   plt.xlabel('Iteration')
   plt.ylabel('Loss')
   plt.title('ICA convergence')
   plt.show()

**Component sign ambiguity:**

.. code-block:: python

   # ICA components have arbitrary sign
   # Ensure consistent orientation
   
   def orient_components(M_matrix, reference_genes):
       """Orient components to have positive weights for reference genes."""
       M_oriented = M_matrix.copy()
       
       for col in M_oriented.columns:
           ref_weights = M_oriented.loc[reference_genes, col].mean()
           if ref_weights < 0:
               M_oriented[col] *= -1
       
       return M_oriented
   
   # Apply to results
   reference_genes = ['housekeeping_gene1', 'housekeeping_gene2']
   M_oriented = {
       species: orient_components(M, reference_genes)
       for species, M in M_matrices.items()
   }