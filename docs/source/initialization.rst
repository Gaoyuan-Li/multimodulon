Initialization of MultiModulon Object
=====================================

This section covers how to initialize and configure the MultiModulon object for your analysis.

Basic Initialization
--------------------

The MultiModulon object is the central class for all analyses:

.. code-block:: python

   from multimodulon import MultiModulon
   
   # Initialize with path to your Input_Data directory
   mm = MultiModulon("path/to/Input_Data")

During initialization, MultiModulon will:

1. Scan the input directory for species/strain subdirectories
2. Load available data files for each species
3. Validate data consistency
4. Set up internal data structures

Constructor Parameters
----------------------

.. py:class:: MultiModulon(input_folder_path)

   Main class for multi-species expression analysis.

   :param str input_folder_path: Path to the Input_Data folder containing species/strain subfolders
   
   :raises ValueError: If the input path doesn't exist or contains no valid species data
   
   **Example:**
   
   .. code-block:: python
      
      # Absolute path
      mm = MultiModulon("/home/user/project/Input_Data")
      
      # Relative path
      mm = MultiModulon("./Input_Data")

Accessing Loaded Data
---------------------

After initialization, you can inspect what was loaded:

View Summary
~~~~~~~~~~~~

.. code-block:: python

   # Print summary of all loaded data
   mm.summary()
   
   # Output:
   # MultiModulon Summary:
   # Species loaded: 3
   # - Species1: 100 samples, 3000 genes
   # - Species2: 80 samples, 2800 genes  
   # - Species3: 120 samples, 3200 genes

List Species
~~~~~~~~~~~~

.. code-block:: python

   # Get list of loaded species
   species_list = mm.species
   print(species_list)
   # ['Species1', 'Species2', 'Species3']

Access Species Data
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Access data for a specific species using dictionary syntax
   species1_data = mm['Species1']
   
   # Or iterate through all species
   for species_name in mm.species:
       data = mm[species_name]
       print(f"{species_name}: {data.log_tpm.shape}")

Species Data Container
----------------------

Each species data is stored in a ``SpeciesData`` object with these attributes:

.. py:class:: SpeciesData

   Container for single species data.
   
   **Attributes (lazy-loaded):**
   
   * **log_tpm** (*pd.DataFrame*) – Log TPM expression matrix
   * **log_tpm_norm** (*pd.DataFrame*) – Normalized log TPM matrix
   * **X** (*pd.DataFrame*) – Aligned expression matrix (defaults to log_tpm_norm)
   * **M** (*pd.DataFrame*) – ICA mixing matrix (gene weights)
   * **A** (*pd.DataFrame*) – ICA activity matrix
   * **sample_sheet** (*pd.DataFrame*) – Sample metadata
   * **gene_table** (*pd.DataFrame*) – Gene annotations
   * **M_thresholds** (*dict*) – Thresholds for M matrix
   * **presence_matrix** (*pd.DataFrame*) – Binarized M matrix

   **Example:**
   
   .. code-block:: python
      
      # Access species data
      species_data = mm['Species1']
      
      # Access expression matrix
      expr_matrix = species_data.log_tpm
      print(f"Expression matrix shape: {expr_matrix.shape}")
      
      # Access sample metadata
      samples = species_data.sample_sheet
      print(f"Conditions: {samples['condition'].unique()}")

Data Validation
---------------

MultiModulon performs automatic validation during initialization:

Manual Validation
~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Validate data for all species
   for species in mm.species:
       is_valid = mm[species].validate_data()
       if not is_valid:
           print(f"Warning: {species} has data inconsistencies")

Common Validation Checks
~~~~~~~~~~~~~~~~~~~~~~~~

1. **Sample consistency**: Sample names in expression matrix match sample sheet
2. **Gene consistency**: Gene IDs are consistent across files  
3. **Data types**: Numeric values in expression matrices
4. **File formats**: Proper CSV formatting

Handling Missing Data
---------------------

MultiModulon gracefully handles missing optional files:

.. code-block:: python

   species_data = mm['Species1']
   
   # Check if optional data is available
   if species_data.gene_table is not None:
       print("Gene annotations available")
   else:
       print("No gene annotations - using gene IDs only")
   
   # Check if ICA has been run
   if species_data.M is not None:
       print("ICA results available")
   else:
       print("ICA not yet performed")

Configuration Options
---------------------

While MultiModulon has sensible defaults, you can configure behavior:

Expression Matrix Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # By default, X uses log_tpm_norm if available, otherwise log_tpm
   # You can manually set which matrix to use for analysis:
   
   species_data = mm['Species1']
   
   # Force use of non-normalized data
   species_data._X = species_data.log_tpm
   
   # Or create custom normalized matrix
   custom_norm = some_normalization_function(species_data.log_tpm)
   species_data._X = custom_norm

Memory Management
~~~~~~~~~~~~~~~~~

For large datasets, data is lazy-loaded to conserve memory:

.. code-block:: python

   # Data is only loaded when first accessed
   species_data = mm['Species1']  # No data loaded yet
   
   # This triggers loading of expression matrix
   expr = species_data.log_tpm  # Now data is loaded and cached
   
   # Subsequent accesses use cached data
   expr2 = species_data.log_tpm  # Returns cached data

Next Steps
----------

After initialization:

1. :doc:`gene_alignment` - Align genes across species
2. :doc:`optimization` - Optimize component numbers  
3. :doc:`multiview_ica` - Run multi-view ICA analysis

Troubleshooting
---------------

**No species loaded:**

.. code-block:: python

   # Check your directory structure
   import os
   print(os.listdir("path/to/Input_Data"))
   # Should show species subdirectories

**Data validation fails:**

.. code-block:: python

   # Check specific issues
   species_data = mm['Species1']
   
   # Verify sample names match
   expr_samples = set(species_data.log_tpm.columns)
   sheet_samples = set(species_data.sample_sheet.index)
   
   missing = expr_samples - sheet_samples
   if missing:
       print(f"Samples in expression but not metadata: {missing}")