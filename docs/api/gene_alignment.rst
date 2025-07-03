Gene Alignment API Reference
============================

This section provides detailed API documentation for gene alignment functions.

BBH Generation
--------------

.. autofunction:: multimodulon.gene_alignment.generate_BBH

   Generate Bidirectional Best Hits (BBH) between all species pairs.

   **Parameters:**
   
   * **multimodulon** (*MultiModulon*) -- MultiModulon instance
   * **output_path** (*str*) -- Directory for BBH results (default: "Output_BBH")
   * **threads** (*int*) -- Number of threads for BLAST (default: 1)
   
   **Process:**
   
   1. Extracts protein sequences from genome files if not available
   2. Runs all-vs-all BLAST for each species pair
   3. Identifies reciprocal best hits
   4. Saves results as CSV files
   
   **Example:**
   
   .. code-block:: python
      
      from multimodulon.gene_alignment import generate_BBH
      
      # Basic usage
      generate_BBH(mm, output_path="Output_BBH", threads=8)
      
      # Check generated files
      import os
      bbh_files = os.listdir("Output_BBH")
      print(f"Generated {len(bbh_files)} BBH files")

Gene Alignment
--------------

.. autofunction:: multimodulon.gene_alignment.align_genes

   Align genes across species using Union-Find algorithm.

   **Parameters:**
   
   * **multimodulon** (*MultiModulon*) -- MultiModulon instance
   * **input_bbh_dir** (*str*) -- BBH directory (default: "Output_BBH")
   * **output_dir** (*str*) -- Output directory (default: "Output_Gene_Info")
   * **reference_order** (*List[str]*) -- Species column order (optional)
   * **bbh_threshold** (*float*) -- Minimum percent identity (optional)
   
   **Returns:**
   
   * **combined_gene_db** (*pd.DataFrame*) -- Aligned gene database
   
   **Example:**
   
   .. code-block:: python
      
      # Default alignment
      gene_db = align_genes(mm)
      
      # With custom parameters
      gene_db = align_genes(
          mm,
          input_bbh_dir="Output_BBH",
          output_dir="Gene_Alignment",
          reference_order=['Species1', 'Species2', 'Species3'],
          bbh_threshold=50  # 50% identity threshold
      )
      
      # Examine results
      print(f"Total gene families: {len(gene_db)}")
      print(f"Core genes: {gene_db.notna().all(axis=1).sum()}")

Union-Find Implementation
-------------------------

.. autoclass:: multimodulon.gene_alignment.UnionFind
   :members:
   :undoc-members:
   
   Union-Find data structure for efficient gene grouping.
   
   **Methods:**
   
   .. automethod:: find
   .. automethod:: union
   
   **Example:**
   
   .. code-block:: python
      
      from multimodulon.gene_alignment import UnionFind
      
      # Create Union-Find structure
      uf = UnionFind()
      
      # Add gene relationships
      uf.union('geneA', 'geneB')  # geneA and geneB are orthologs
      uf.union('geneB', 'geneC')  # geneB and geneC are orthologs
      
      # Find gene families
      print(uf.find('geneA'))  # Returns root of family
      print(uf.find('geneC'))  # Same root as geneA

BBH Utilities
-------------

BBHAnalyzer Class
~~~~~~~~~~~~~~~~~

.. autoclass:: multimodulon.utils.bbh.BBHAnalyzer
   :members:
   :undoc-members:
   
   Performs Bidirectional Best Hits analysis using containerized BLAST.
   
   **Constructor:**
   
   .. automethod:: __init__
   
   **Methods:**
   
   .. automethod:: run_bbh_all_pairs
   .. automethod:: run_blastp
   .. automethod:: parse_blast_results
   .. automethod:: find_bbh_pairs
   
   **Example:**
   
   .. code-block:: python
      
      from multimodulon.utils.bbh import BBHAnalyzer
      
      # Initialize with custom Docker image
      analyzer = BBHAnalyzer(
          docker_image="quay.io/biocontainers/blast:2.16.0--h66d330f_5"
      )
      
      # Run BBH for species pairs
      species_data = {
          'Species1': {'protein_path': 'path/to/species1.faa'},
          'Species2': {'protein_path': 'path/to/species2.faa'}
      }
      
      bbh_results = analyzer.run_bbh_all_pairs(species_data)
      
      # Access results
      for (sp1, sp2), bbh_df in bbh_results.items():
          print(f"{sp1} vs {sp2}: {len(bbh_df)} BBH pairs")

Helper Functions
----------------

Creating Gene Tables
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: multimodulon.gene_alignment.create_gene_tables

   Create gene annotation tables from GFF files.

   **Parameters:**
   
   * **species_paths** (*Dict[str, Path]*) -- Paths to species directories
   * **gff_feature** (*str*) -- Feature type to extract (default: "CDS")
   
   **Example:**
   
   .. code-block:: python
      
      species_paths = {
          'Species1': Path('Input_Data/Species1'),
          'Species2': Path('Input_Data/Species2')
      }
      
      create_gene_tables(species_paths, gff_feature="CDS")

Processing BBH Results
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: multimodulon.gene_alignment.process_bbh_files

   Process multiple BBH result files.

   **Parameters:**
   
   * **bbh_dir** (*str*) -- Directory containing BBH CSV files
   * **threshold** (*float*) -- Minimum percent identity (optional)
   
   **Returns:**
   
   * **all_pairs** (*List[Tuple[str, str]]*) -- List of ortholog pairs
   
   **Example:**
   
   .. code-block:: python
      
      # Get all ortholog pairs above threshold
      ortholog_pairs = process_bbh_files(
          bbh_dir="Output_BBH",
          threshold=60  # 60% identity
      )
      
      print(f"Total ortholog pairs: {len(ortholog_pairs)}")

Usage Examples
--------------

Complete Gene Alignment Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from multimodulon import MultiModulon
   from multimodulon.gene_alignment import generate_BBH, align_genes
   
   # Initialize
   mm = MultiModulon("Input_Data")
   
   # Step 1: Generate BBH
   print("Generating BBH...")
   generate_BBH(mm, output_path="BBH_Results", threads=16)
   
   # Step 2: Align genes with different thresholds
   thresholds = [30, 40, 50, 60, 70]
   alignment_stats = {}
   
   for threshold in thresholds:
       print(f"\nAligning with {threshold}% identity threshold...")
       gene_db = align_genes(
           mm,
           input_bbh_dir="BBH_Results",
           output_dir=f"Alignment_{threshold}",
           bbh_threshold=threshold
       )
       
       # Calculate statistics
       total_families = len(gene_db)
       core_genes = gene_db.notna().all(axis=1).sum()
       alignment_stats[threshold] = {
           'total': total_families,
           'core': core_genes,
           'coverage': gene_db.notna().sum().sum() / gene_db.size
       }
   
   # Display results
   import pandas as pd
   stats_df = pd.DataFrame(alignment_stats).T
   print("\nAlignment Statistics:")
   print(stats_df)

Custom BBH Analysis
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from multimodulon.utils.bbh import BBHAnalyzer
   import subprocess
   
   # Custom BLAST parameters
   class CustomBBHAnalyzer(BBHAnalyzer):
       def run_blastp(self, query_fasta, subject_fasta, output_file, threads=1):
           """Run BLASTP with custom parameters."""
           cmd = [
               'blastp',
               '-query', query_fasta,
               '-subject', subject_fasta,
               '-out', output_file,
               '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
               '-num_threads', str(threads),
               '-max_target_seqs', '1',  # Only top hit
               '-evalue', '1e-10',       # Stricter e-value
               '-qcov_hsp_perc', '80'    # Query coverage threshold
           ]
           subprocess.run(cmd, check=True)
   
   # Use custom analyzer
   analyzer = CustomBBHAnalyzer()
   bbh_results = analyzer.run_bbh_all_pairs(species_data)

Incremental Gene Alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Add new species to existing alignment
   existing_gene_db = pd.read_csv("Output_Gene_Info/combined_gene_db.csv", index_col=0)
   existing_species = list(existing_gene_db.columns)
   
   # Generate BBH only for new species pairs
   new_species = 'SpeciesNew'
   for species in existing_species:
       print(f"Running BBH: {new_species} vs {species}")
       # Run BBH for this pair
       # Save to BBH directory
   
   # Re-run alignment with all BBH files
   updated_gene_db = align_genes(
       mm,
       reference_order=existing_species + [new_species]
   )

Analyzing Ortholog Relationships
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load gene database
   gene_db = pd.read_csv("Output_Gene_Info/combined_gene_db.csv", index_col=0)
   
   # Find multi-copy orthologs
   for species in gene_db.columns:
       gene_counts = gene_db[species].value_counts()
       multi_copy = gene_counts[gene_counts > 1]
       if len(multi_copy) > 0:
           print(f"\n{species} multi-copy genes:")
           print(multi_copy.head())
   
   # Create ortholog network
   import networkx as nx
   
   G = nx.Graph()
   
   # Add edges for orthologs
   for idx, row in gene_db.iterrows():
       genes = row.dropna().tolist()
       for i in range(len(genes)):
           for j in range(i+1, len(genes)):
               G.add_edge(genes[i], genes[j])
   
   # Analyze network
   print(f"Ortholog network: {G.number_of_nodes()} genes, {G.number_of_edges()} relationships")
   
   # Find connected components
   components = list(nx.connected_components(G))
   print(f"Number of gene families: {len(components)}")
   
   # Largest families
   largest = sorted(components, key=len, reverse=True)[:5]
   for i, family in enumerate(largest):
       print(f"Family {i+1}: {len(family)} genes")