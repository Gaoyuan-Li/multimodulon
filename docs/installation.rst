Installation
============

This guide will help you install MultiModulon and its dependencies.

Requirements
------------

* Python 3.8 or higher
* CUDA-capable GPU (optional, for GPU acceleration)
* Docker (optional, for BBH analysis using containerized BLAST)

Basic Installation
------------------

Install MultiModulon using pip:

.. code-block:: bash

   pip install multimodulon

Development Installation
------------------------

For development or to get the latest features:

.. code-block:: bash

   git clone https://github.com/yourusername/multimodulon.git
   cd multimodulon
   pip install -e .

Dependencies
------------

Core dependencies will be automatically installed:

* numpy
* pandas
* matplotlib
* scikit-learn
* pytorch (for multi-view ICA)
* biopython (for sequence handling)

Optional Dependencies
---------------------

For GPU support:

.. code-block:: bash

   pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

For BBH analysis with Docker:

.. code-block:: bash

   docker pull quay.io/biocontainers/blast:2.16.0--h66d330f_5

Verifying Installation
----------------------

To verify that MultiModulon is installed correctly:

.. code-block:: python

   import multimodulon
   print(multimodulon.__version__)

Troubleshooting
---------------

If you encounter issues:

1. **ImportError**: Make sure all dependencies are installed
2. **CUDA errors**: Check that your PyTorch installation matches your CUDA version
3. **Docker errors**: Ensure Docker daemon is running and you have proper permissions

For more help, please open an issue on our GitHub repository.