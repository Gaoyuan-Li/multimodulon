.. MultiModulon documentation master file

Welcome to MultiModulon's documentation!
========================================

**MultiModulon** is a Python package for an integrative multi-species/multi-strains/multi-modalities framework. It enables the identification of conserved and species-specific regulatory modules across multiple bacterial strains or species.

.. note::
   This project is under active development. The API may change in future releases.

Key Features
------------

* **Multi-species Integration**: Analyze gene expression data from multiple bacterial species simultaneously
* **Multi-view ICA**: Identify core (conserved) and unique (species-specific) regulatory modules
* **Optimization Tools**: Automated optimization of component numbers using multiple metrics
* **Visualization**: Comprehensive plotting functions for iModulon activities and gene weights
* **GPU Support**: Accelerated computation using GPU for large-scale analyses

.. note::
   This project is under active development. The API may change in future releases.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
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