Installation
============

This guide will help you install MultiModulon and its dependencies.

Requirements
------------

* Python 3.10 or higher
* CUDA-capable GPU
* Docker

Installation
------------------

Install BLAST+ (Required for BBH generation)

.. code-block:: bash

   conda install -c bioconda blast

Install PyTorch (Required for multi-view ICA)

.. code-block:: bash

   # Install PyTorch with CUDA support
   pip install torch==2.6.0 torchvision==0.21.0 torchaudio==2.6.0 --index-url https://download.pytorch.org/whl/cu124

   # Install geotorch
   pip install geotorch==0.3.0


Install MultiModulon using pip:

.. code-block:: bash

   pip install git+https://github.com/Gaoyuan-Li/multimodulon.git

Development Installation
------------------------

For development:

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
* pytorch (for robust multi-view ICA)
* biopython (for sequence handling)


Verifying Installation
----------------------

To verify that MultiModulon is installed correctly:

.. code-block:: python

   import multimodulon
   print(multimodulon.__version__)
