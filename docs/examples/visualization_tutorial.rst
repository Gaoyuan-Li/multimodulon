Characterization of Core Components
===================================

This notebook demonstrates the second step for multi-species/strain/modality analysis using the MultiModulon package.

Step 1: Load MultiModulon object
--------------------------------

Load the saved MultiModulon object from Tutorial_1_Create_multiModulon_object.ipynb.

.. code-block:: python

   # Import required libraries
   from multimodulon import MultiModulon
   import pandas as pd
   import numpy as np
   import matplotlib.pyplot as plt

   # Set display options
   pd.set_option('display.max_columns', None)
   pd.set_option('display.max_rows', 50)

   multiModulon = MultiModulon.load_json_multimodulon("./multiModulon_E_coli_comparison_demo.json.gz")

   multiModulon.species_palette

Step 2: Visualize the Core Components (Available Functions)
------------------------------------------------------------

Visualize the core and unique components by ``view_iModulon_weights`` and ``view_core_iModulon_weights``.

.. code-block:: python

   multiModulon.view_iModulon_weights(
       species="MG1655",
       component="Core_1",
       fig_size=(8, 5),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       show_COG=True,
       show_gene_names=True
   )

.. code-block:: python

   multiModulon.view_iModulon_genes(
       species="MG1655",
       component="Core_1"
   )

.. code-block:: python

   multiModulon['MG1655'].presence_matrix

.. code-block:: python

   multiModulon['MG1655'].sample_sheet

.. code-block:: python

   multiModulon.view_iModulon_activities(
       species="MG1655",
       component="Core_1",
       fig_size=(6, 4),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
   )

.. code-block:: python

   multiModulon.compare_core_iModulon_activity(
       component='Core_1',
       fig_size=(4, 3),
       species_in_comparison=['MG1655', 'BL21', 'W3110'],
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       condition_list=['M9_glucose:in_house_2', 'LiAcet:in_house_2', 'acetate:in_house_2'],
       legend_title='Strain',
       title=None
   )

.. code-block:: python

   multiModulon.show_iModulon_activity_change(
       species='MG1655',
       condition_1='acetate:in_house_2',
       condition_2='LiAcet:in_house_2',
       fig_size=(5, 5),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
   )

.. code-block:: python

   multiModulon.show_gene_iModulon_correlation(
       gene='b0583',
       component='Core_1',
       fig_size=(4, 3),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
   )

Step 3: Compare Gene Membership for Core Components
----------------------------------------------------

Visualize gene weights and compare gene membership across strains.

.. code-block:: python

   multiModulon.view_core_iModulon_weights(
       component="Core_1",
       fig_size=(6, 4),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       reference_order=['MG1655', 'BL21', 'W3110'],
       show_COG=True
   )

.. code-block:: python

   Core_1_comparison = multiModulon.compare_core_iModulon(
       component='Core_1',
       y_label='Strains',
       reference_order=['MG1655', 'BL21', 'W3110'],
       fig_size=(20, 6),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       show_gene_names=False
       # save_path="./Output_iModulon_Figures/"
   )

.. code-block:: python

   Core_1_comparison
