Optimization of Thresholds
==========================

This section covers the optimization of thresholds for iModulon gene weights using Otsu's method, which automatically determines optimal cutoffs to distinguish significant genes from background noise.

Overview
--------

After running multi-view ICA, gene weights in the M matrices need to be thresholded to identify which genes significantly contribute to each iModulon. The ``optimize_M_thresholds`` method applies an adapted version of Otsu's method to automatically determine these thresholds.

What is Otsu's Method?
----------------------

`Otsu's method <https://en.wikipedia.org/wiki/Otsu%27s_method>`_ is an automatic threshold selection technique originally developed for image segmentation. The method finds the optimal threshold by maximizing the between-class variance, effectively separating data into two distinct groups.

**Key principles:**

* Analyzes the distribution of values (histogram)
* Divides data into two classes: signal and background
* Finds the threshold that maximizes separation between classes
* Works best when data has a bimodal distribution

**Mathematical basis:**

The method maximizes the between-class variance:

.. math::

   \sigma_b^2(t) = \omega_0(t)\omega_1(t)[\mu_0(t) - \mu_1(t)]^2

Where:
   - :math:`\omega_0` and :math:`\omega_1` are the class probabilities
   - :math:`\mu_0` and :math:`\mu_1` are the class mean values
   - :math:`t` is the threshold being evaluated

Adaptation for Gene Expression Data
------------------------------------

Gene expression data often has heavy-tailed distributions with many genes showing low weights and few showing high weights. To handle this, MultiModulon uses a **quantile-based pre-filtering approach**:

1. **Pre-filtering**: Remove the bottom X% of absolute values (default: 90%)
2. **Apply Otsu**: Run Otsu's method on the remaining high-value genes
3. **Binarization**: Create presence matrix based on optimized thresholds

This approach helps identify true regulatory genes while filtering out noise.

Basic Usage
-----------

.. py:method:: MultiModulon.optimize_M_thresholds(method="Otsu's method", quantile_threshold=90)

   Optimize thresholds for M matrices across all species.

   :param str method: Threshold optimization method (currently only "Otsu's method")
   :param float quantile_threshold: Percentile for pre-filtering (default: 90)

**Example:**

.. code-block:: python

   # After running ICA
   multiModulon.run_robust_multiview_ica(
       a={'Species1': 50, 'Species2': 60},
       c=20,
       num_runs=100
   )
   
   # Optimize thresholds
   multiModulon.optimize_M_thresholds(
       method="Otsu's method",
       quantile_threshold=90
   )
   
   # Results are stored in each species
   for species in multiModulon.species:
       thresholds = multiModulon[species].M_thresholds
       presence = multiModulon[species].presence_matrix
       print(f"{species}: {presence.sum().mean():.1f} genes per component")

Understanding the Results
-------------------------

M_thresholds DataFrame
~~~~~~~~~~~~~~~~~~~~~~

Contains the optimized threshold for each component:

.. code-block:: python

   # Access thresholds
   thresholds = multiModulon['Species1'].M_thresholds
   
   # Structure:
   # Index: Component names (Core_1, Core_2, ..., Unique_1, ...)
   # Column: 'M_threshold' - the optimized threshold value
   
   print("Component thresholds:")
   print(thresholds.head(10))
   
   # Average threshold across components
   avg_threshold = thresholds['M_threshold'].mean()
   print(f"Average threshold: {avg_threshold:.3f}")

Presence Matrix
~~~~~~~~~~~~~~~

Binary matrix indicating which genes belong to each component:

.. code-block:: python

   # Access presence matrix
   presence = multiModulon['Species1'].presence_matrix
   
   # Structure:
   # Rows: Genes
   # Columns: Components
   # Values: 1 if |gene weight| > threshold, 0 otherwise
   
   # Count genes per component
   genes_per_comp = presence.sum(axis=0)
   print(f"Genes per component: min={genes_per_comp.min()}, max={genes_per_comp.max()}")
   
   # Find components with many genes
   large_components = genes_per_comp[genes_per_comp > 100].index
   print(f"Components with >100 genes: {large_components.tolist()}")
   
   # Get genes in a specific component
   component = 'Core_1'
   component_genes = presence[presence[component] == 1].index.tolist()
   print(f"{component} contains {len(component_genes)} genes")

Customizing Pre-filtering
-------------------------

Adjust the quantile threshold based on your data:

.. code-block:: python

   # More stringent: Keep only top 5% of values
   multiModulon.optimize_M_thresholds(quantile_threshold=95)
   
   # Less stringent: Keep top 20% of values  
   multiModulon.optimize_M_thresholds(quantile_threshold=80)
   
   # Compare results
   for threshold in [80, 85, 90, 95]:
       multiModulon.optimize_M_thresholds(quantile_threshold=threshold)
       avg_genes = []
       for species in multiModulon.species:
           presence = multiModulon[species].presence_matrix
           avg_genes.append(presence.sum().mean())
       print(f"Quantile {threshold}%: {np.mean(avg_genes):.1f} genes/component")

Next Steps
----------

After optimizing thresholds:

1. :doc:`visualization` - Visualize thresholded components
2. Use presence matrices for enrichment analysis
3. Compare gene sets across species for core components
4. Export gene lists for biological validation