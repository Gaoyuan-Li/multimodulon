.. MultiModulon documentation master file

Welcome to MultiModulon's documentation!
========================================

**MultiModulon** is a Python package for an integrative multi-species/multi-strains/multi-modalities analysis framework. It enables the identification of conserved and species-specific regulatory modules across multiple bacterial strains or species.

.. image:: https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54
   :target: https://www.python.org
   :alt: Python 3.10+

.. image:: https://img.shields.io/badge/PyTorch-EE4C2C?style=for-the-badge&logo=pytorch&logoColor=white
   :target: https://pytorch.org
   :alt: PyTorch 2.6+

.. image:: https://img.shields.io/badge/Matplotlib-%23ffffff.svg?style=for-the-badge&logo=Matplotlib&logoColor=black
   :target: https://matplotlib.org/stable/users/index.html
   :alt: Matplotlib

.. image:: https://img.shields.io/badge/pandas-%23150458.svg?style=for-the-badge&logo=pandas&logoColor=white
   :target: https://pandas.pydata.org/
   :alt: Pandas

.. image:: https://img.shields.io/badge/scikit--learn-%23F7931E.svg?style=for-the-badge&logo=scikit-learn&logoColor=white
   :target: https://scikit-learn.org/stable/
   :alt: scikit-learn

.. image:: https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white
   :target: https://www.docker.com/
   :alt: Docker

Key Features
------------

* **Multi-modalities/Multi-species Integration**: Analyze gene expression data from multiple bacterial species simultaneously
* **Robust Multi-view ICA**: Identify core (conserved) and unique (species-specific) regulatory modules
* **Optimization Tools**: Automated optimization of component numbers using non-single-gene filtering
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
   :caption: Examples & Tutorials
   
   examples/basic_workflow
   examples/characterization_core_components
   examples/characterization_unique_components

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
   :maxdepth: 1
   :caption: Additional Information
   
   changelog
   error_report

.. note::
   This project is under active development. The API may change in future releases.
