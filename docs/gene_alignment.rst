Gene Alignment
==============

This section covers gene alignment across multiple species using Bidirectional Best Hits (BBH) and the creation of a unified gene database.

Overview
--------

Gene alignment is a crucial step that enables multi-species analysis by:

1. Identifying orthologous genes across species
2. Creating a unified gene indexing system
3. Aligning expression matrices for multi-view ICA

The process involves two main steps:

1. **BBH Generation**: Find ortholog pairs between all species pairs
2. **Gene Alignment**: Group orthologs into gene families using Union-Find

Generating BBH
--------------

Basic BBH Generation
~~~~~~~~~~~~~~~~~~~~

.. py:method:: MultiModulon.generate_BBH(output_path="Output_BBH", threads=1)

   Generate Bidirectional Best Hits between all species pairs.

   :param str output_path: Directory to save BBH results (default: "Output_BBH")
   :param int threads: Number of threads for BLAST (default: 1)
   
   **Example:**
   
   .. code-block:: python
      
      # Generate BBH with multiple threads
      multiModulon.generate_BBH(
          output_path="Output_BBH",
          threads=8
      )
      
      # Check results
      print(f"BBH files generated in Output_BBH/")

How BBH Works
~~~~~~~~~~~~~

1. **Protein sequences**: Extracts from genome files or uses provided protein.faa
2. **BLAST search**: All-vs-all BLAST between each species pair
3. **Best hits**: Identifies reciprocal best hits
4. **Output format**: CSV files with columns: query, subject, pident, length, etc.

Requirements
~~~~~~~~~~~~

For BBH generation, each species needs either:

* **Option 1**: ``genome.fna`` and ``genome.gff`` files
* **Option 2**: Pre-computed ``protein.faa`` file

Docker Support
~~~~~~~~~~~~~~

BBH uses containerized BLAST for reproducibility:

.. code-block:: python

   # The default Docker image is used automatically
   # To use a different BLAST version:
   from multimodulon.utils.bbh import BBHAnalyzer
   
   analyzer = BBHAnalyzer(
       docker_image="quay.io/biocontainers/blast:2.16.0--h66d330f_5"
   )

Aligning Genes
--------------

After BBH generation, align genes across species:

.. py:method:: MultiModulon.align_genes(input_bbh_dir="Output_BBH", output_dir="Output_Gene_Info", reference_order=None, bbh_threshold=None)

   Align genes across all species using Union-Find algorithm.

   :param str input_bbh_dir: Directory containing BBH CSV files
   :param str output_dir: Directory to save combined gene database  
   :param list reference_order: Order of species columns in output (optional)
   :param float bbh_threshold: Minimum percent identity for BBH (optional)
   :return: Combined gene database DataFrame
   :rtype: pd.DataFrame
   
   **Example:**
   
   .. code-block:: python
      
      # Basic alignment
      combined_gene_db = multiModulon.align_genes() 
      
      # With custom parameters
      combined_gene_db = multiModulon.align_genes(
          input_bbh_dir="Output_BBH",
          output_dir="Output_Gene_Info",
          reference_order=['Species1', 'Species2', 'Species3'],
          bbh_threshold=90  # 90% identity threshold
      )
      
      # Examine the results
      print(combined_gene_db.head())

Combined Gene Database Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The output is a DataFrame where:

* **Rows**: Gene families (groups of orthologs)
* **Columns**: Species names
* **Values**: Gene IDs for each species (NaN if absent)

Example output:

.. code-block:: text

   Species1    Species2    Species3    row_label
   gene001     geneA_001   locus_001   gene001
   gene002     geneA_002   NaN         gene002
   NaN         geneA_003   locus_003   geneA_003
   NaN         NaN         locus_004   locus_004

Union-Find Algorithm
~~~~~~~~~~~~~~~~~~~~

The alignment uses Union-Find to group genes:

1. Each gene starts in its own group
2. BBH relationships merge groups
3. Transitive closure creates gene families
4. Result: Genes in same family are orthologs

Filtering Options
~~~~~~~~~~~~~~~~~

Control alignment stringency:

.. code-block:: python

   # Strict alignment - high identity threshold
   strict_db = multiModulon.align_genes(bbh_threshold=80)
   
   # Permissive alignment - lower threshold  
   permissive_db = multiModulon.align_genes(bbh_threshold=30)
   
   # Check alignment statistics
   print(f"Strict: {strict_db.notna().sum().sum()} genes aligned")
   print(f"Permissive: {permissive_db.notna().sum().sum()} genes aligned")

Generating Expression Matrices
------------------------------

After alignment, generate aligned expression matrices:

.. py:method:: MultiModulon.generate_X(gene_info_folder)

   Generate X matrices with consistent row indices based on combined_gene_db.

   :param str gene_info_folder: Path to folder containing combined_gene_db.csv
   
   **Example:**
   
   .. code-block:: python
      
      # Generate aligned expression matrices
      multiModulon.generate_X("Output_Gene_Info")
      
      # Access aligned matrices
      for species in multiModulon.species:
          X = multiModulon[species].X
          print(f"{species}: {X.shape}")
          print(f"Missing genes: {X.isna().sum().sum()}")

The aligned matrices have:

* **Consistent row order**: Same gene families in same positions
* **Ready for ICA**: Can be directly used for multi-view ICA

Quality Control
---------------

Check alignment quality:

.. code-block:: python

   # Alignment statistics
   total_families = len(combined_gene_db)
   
   # Genes per species
   for species in combined_gene_db.columns:
       gene_count = combined_gene_db[species].notna().sum()
       print(f"{species}: {gene_count} genes")
   
   # Core genes (present in all species)
   core_genes = combined_gene_db.notna().all(axis=1).sum()
   print(f"Core genes: {core_genes}")
   
   # Species-specific genes
   for species in combined_gene_db.columns:
       specific = (
           combined_gene_db[species].notna() & 
           combined_gene_db.drop(columns=species).isna().all(axis=1)
       ).sum()
       print(f"{species}-specific: {specific}")

Advanced Usage
--------------

Custom Gene Grouping
~~~~~~~~~~~~~~~~~~~~

For custom ortholog definitions:

.. code-block:: python

   # Load your own ortholog mappings
   custom_orthologs = pd.read_csv("custom_orthologs.csv")
   
   # Create combined gene database manually
   from multimodulon.gene_alignment import create_combined_gene_db
   combined_db = create_combined_gene_db(
       custom_orthologs,
       species_list=multiModulon.species   
   )

Next Steps
----------

After gene alignment:

1. :doc:`optimization` - Optimize component numbers
2. :doc:`multiview_ica` - Run multi-view ICA
3. :doc:`visualization` - Visualize aligned components