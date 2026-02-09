Characterization of Unique Components
=====================================

This notebook demonstrates the third step for multi-species/strain/modality analysis using the MultiModulon package.

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

Step 2: Characterize Unique Components
--------------------------------------

Characterize unique components for each species.

MG1655
~~~~~~

.. code-block:: python

   multiModulon.view_iModulon_weights(
       species="MG1655",
       component="Unique_1",
       fig_size=(6, 4),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       show_COG=True
   )

.. code-block:: python

   multiModulon.view_iModulon_activities(
       species="MG1655",
       component="Unique_1",
       fig_size=(8, 3),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
   )

BL21
~~~~

.. code-block:: python

   multiModulon.view_iModulon_weights(
       species="BL21",
       component="Unique_1",
       fig_size=(6, 4),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       show_COG=True
   )

.. code-block:: python

   multiModulon.view_iModulon_activities(
       species="BL21",
       component="Unique_1",
       fig_size=(8, 3),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       highlight_project=['in_house_1', 'multi_ALE', 'pan_fur']
   )

.. code-block:: python

   multiModulon.view_iModulon_weights(
       species="BL21",
       component="Unique_3",
       fig_size=(6, 4),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       show_COG=True
   )

.. code-block:: python

   multiModulon.view_iModulon_activities(
       species="BL21",
       component="Unique_3",
       fig_size=(8, 3),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       highlight_project=['in_house_2']
   )

W3110
~~~~~

.. code-block:: python

   multiModulon.view_iModulon_weights(
       species="W3110",
       component="Unique_2",
       fig_size=(6, 4),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
       show_COG=True
   )

.. code-block:: python

   multiModulon.view_iModulon_activities(
       species="W3110",
       component="Unique_2",
       fig_size=(8, 3),
       font_path="/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
   )
