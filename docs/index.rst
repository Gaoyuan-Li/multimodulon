.. MultiModulon documentation master file

Welcome to MultiModulon's documentation!
========================================

**MultiModulon** is a Python package for integrative multi-species expression analysis using multi-view Independent Component Analysis (ICA). It enables the identification of conserved and species-specific regulatory modules across multiple bacterial strains or species.

Key Features
------------

* **Multi-species Integration**: Analyze gene expression data from multiple bacterial species simultaneously
* **Gene Alignment**: Automatic ortholog detection and gene alignment across species using Bidirectional Best Hits (BBH)
* **Multi-view ICA**: Identify core (conserved) and unique (species-specific) regulatory modules
* **Optimization Tools**: Automated optimization of component numbers using multiple metrics
* **Visualization**: Comprehensive plotting functions for iModulon activities and gene weights
* **GPU Support**: Accelerated computation using GPU for large-scale analyses

.. note::
   This project is under active development. The API may change in future releases.

Quick Start
-----------

.. code-block:: python

   from multimodulon import MultiModulon
   
   # Initialize with your data directory
   mm = MultiModulon("path/to/Input_Data")
   
   # Generate BBH and align genes
   mm.generate_BBH(threads=8)
   mm.align_genes()
   
   # Run multi-view ICA
   mm.run_multiview_ica(a={'species1': 50, 'species2': 60}, c=20)
   
   # Visualize results
   mm.view_iModulon_activities('species1', 'Core_1')

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   data_preparation
   
.. toctree::
   :maxdepth: 2
   :caption: Core Functionality
   
   initialization
   gene_alignment
   optimization
   multiview_ica
   visualization
   
.. toctree::
   :maxdepth: 2
   :caption: API Reference
   
   api/core
   api/optimization
   api/gene_alignment
   api/multiview_ica
   api/visualization
   api/utilities
   
.. toctree::
   :maxdepth: 2
   :caption: Examples & Tutorials
   
   examples/basic_workflow
   examples/optimization_tutorial
   examples/visualization_tutorial
   
.. toctree::
   :maxdepth: 1
   :caption: Additional Information
   
   changelog
   contributing