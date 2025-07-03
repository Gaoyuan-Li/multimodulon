.. MultiModulon documentation master file

Welcome to MultiModulon's documentation!
========================================

**MultiModulon** is a Python package for an integrative multi-species/multi-strains/multi-modalities analysis framework. It enables the identification of conserved and species-specific regulatory modules across multiple bacterial strains or species.

.. image:: https://img.shields.io/pypi/pyversions/multimodulon?style=for-the-badge&color=85A2C1
   :target: https://www.python.org
   :alt: Python 3.10+

.. image:: https://img.shields.io/badge/PyTorch-EE4C2C?style=for-the-badge&logo=pytorch&logoColor=white
   :target: https://pytorch.org
   :alt: PyTorch 2.6+

.. note::
   This project is under active development. The API may change in future releases.

Key Features
------------

* **Multi-species Integration**: Analyze gene expression data from multiple bacterial species simultaneously
* **Multi-view ICA**: Identify core (conserved) and unique (species-specific) regulatory modules
* **Optimization Tools**: Automated optimization of component numbers using multiple metrics
* **Visualization**: Comprehensive plotting functions for iModulon activities and gene weights
* **GPU Support**: Accelerated computation using GPU for large-scale analyses

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
   threshold_optimization
   visualization
   
.. toctree::
   :maxdepth: 2
   :caption: Examples & Tutorials
   
   examples/basic_workflow
   examples/visualization_tutorial
   
.. toctree::
   :maxdepth: 1
   :caption: Additional Information
   
   changelog
   error_report