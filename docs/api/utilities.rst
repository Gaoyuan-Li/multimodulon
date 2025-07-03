Utilities API Reference
=======================

This section provides API documentation for utility functions and helper modules.

GFF Utilities
-------------

.. autofunction:: multimodulon.gff_utils.gff2pandas

   Convert GFF file(s) to pandas DataFrame with enhanced attribute parsing.

   **Parameters:**
   
   * **gff_file** (*str or list*) -- Path(s) to GFF file(s)
   * **feature** (*str or list*) -- Feature type(s) to extract (default: "CDS")
   * **index** (*str, optional*) -- Column or attribute to use as index
   
   **Returns:**
   
   * **gff_df** (*pd.DataFrame*) -- GFF data as DataFrame
   
   **Features:**
   
   * Handles multiple GFF files
   * Parses all GFF attributes into separate columns
   * Filters by feature type
   * Custom indexing options
   
   **Example:**
   
   .. code-block:: python
      
      from multimodulon.gff_utils import gff2pandas
      
      # Single file
      gff_df = gff2pandas('genome.gff', feature='CDS', index='locus_tag')
      
      # Multiple files
      gff_df = gff2pandas(
          ['genome1.gff', 'genome2.gff'],
          feature=['CDS', 'tRNA', 'rRNA']
      )
      
      # Access parsed attributes
      print(gff_df.columns)  # Shows all parsed attributes

.. autofunction:: multimodulon.gff_utils.create_gene_table

   Create gene annotation tables for all species from GFF files.

   **Parameters:**
   
   * **multimodulon** (*MultiModulon*) -- MultiModulon instance
   
   **Process:**
   
   1. Reads GFF files from each species directory
   2. Extracts CDS features
   3. Creates gene_table.csv with annotations
   
   **Example:**
   
   .. code-block:: python
      
      from multimodulon.gff_utils import create_gene_table
      
      # Create gene tables for all species
      create_gene_table(mm)
      
      # Access created tables
      gene_table = mm['Species1'].gene_table
      print(gene_table.head())

FASTA Utilities
---------------

.. autofunction:: multimodulon.utils.fasta_utils.extract_protein_sequences

   Extract protein sequences from genome using GFF annotations.

   **Parameters:**
   
   * **genome_fasta** (*Path*) -- Path to genome FASTA file
   * **gff_file** (*Path*) -- Path to GFF annotation file
   * **output_fasta** (*Path*) -- Path for output protein FASTA
   
   **Returns:**
   
   * **output_path** (*Path*) -- Path to created protein file
   
   **Features:**
   
   * Extracts CDS sequences
   * Translates to protein
   * Handles multiple genetic codes
   * Adds proper FASTA headers
   
   **Example:**
   
   .. code-block:: python
      
      from pathlib import Path
      from multimodulon.utils.fasta_utils import extract_protein_sequences
      
      # Extract proteins
      protein_file = extract_protein_sequences(
          genome_fasta=Path('genome.fna'),
          gff_file=Path('genome.gff'),
          output_fasta=Path('proteins.faa')
      )
      
      print(f"Proteins saved to: {protein_file}")

BBH Utilities
-------------

.. autoclass:: multimodulon.utils.bbh.BBHAnalyzer

   Bidirectional Best Hits analyzer using containerized BLAST.

   **Constructor:**
   
   .. code-block:: python
      
      BBHAnalyzer(docker_image="quay.io/biocontainers/blast:2.16.0--h66d330f_5")
   
   **Parameters:**
   
   * **docker_image** (*str*) -- Docker image with BLAST tools

.. automethod:: multimodulon.utils.bbh.BBHAnalyzer.run_bbh_all_pairs

   Run BBH analysis for all species pairs.

   **Parameters:**
   
   * **species_data** (*Dict*) -- Species names and protein paths
   
   **Returns:**
   
   * **bbh_results** (*Dict*) -- BBH DataFrames for each pair
   
   **Example:**
   
   .. code-block:: python
      
      from multimodulon.utils.bbh import BBHAnalyzer
      
      analyzer = BBHAnalyzer()
      
      species_data = {
          'Species1': {'protein_path': 'species1.faa'},
          'Species2': {'protein_path': 'species2.faa'},
          'Species3': {'protein_path': 'species3.faa'}
      }
      
      bbh_results = analyzer.run_bbh_all_pairs(species_data)
      
      # Access results
      for (sp1, sp2), bbh_df in bbh_results.items():
          print(f"{sp1} vs {sp2}: {len(bbh_df)} BBH pairs")
          print(f"Average identity: {bbh_df['pident'].mean():.1f}%")

GFF Parser
----------

.. autofunction:: multimodulon.utils.gff_parser.parse_gff

   Parse GFF file and extract gene information.

   **Parameters:**
   
   * **gff_file** (*Path*) -- Path to GFF file
   * **feature_type** (*str*) -- Feature to extract (default: "CDS")
   
   **Returns:**
   
   * **gene_df** (*pd.DataFrame*) -- Gene information indexed by locus_tag
   
   **Extracted Fields:**
   
   * locus_tag
   * gene_name  
   * product
   * start, end, strand
   * Additional GFF attributes
   
   **Example:**
   
   .. code-block:: python
      
      from pathlib import Path
      from multimodulon.utils.gff_parser import parse_gff
      
      # Parse GFF
      genes = parse_gff(Path('genome.gff'))
      
      # Access gene info
      print(f"Total genes: {len(genes)}")
      print(genes[['gene_name', 'product', 'start', 'end']].head())

Expression Data Utilities
-------------------------

Normalization Functions
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def normalize_expression(expr_matrix, method='quantile'):
       """
       Normalize expression matrix.
       
       Parameters
       ----------
       expr_matrix : pd.DataFrame
           Expression matrix (genes x samples)
       method : str
           Normalization method ('quantile', 'tmm', 'median')
       
       Returns
       -------
       pd.DataFrame
           Normalized expression matrix
       """
       if method == 'quantile':
           # Quantile normalization
           from sklearn.preprocessing import QuantileTransformer
           qt = QuantileTransformer(output_distribution='normal')
           normalized = qt.fit_transform(expr_matrix.T).T
           return pd.DataFrame(normalized, 
                             index=expr_matrix.index,
                             columns=expr_matrix.columns)
       
       elif method == 'median':
           # Median normalization
           medians = expr_matrix.median(axis=0)
           target_median = medians.median()
           factors = target_median / medians
           return expr_matrix * factors
       
       else:
           raise ValueError(f"Unknown method: {method}")

Quality Control
~~~~~~~~~~~~~~~

.. code-block:: python

   def qc_expression_matrix(expr_matrix, min_expression=1, min_samples=10):
       """
       Quality control for expression matrix.
       
       Parameters
       ----------
       expr_matrix : pd.DataFrame
           Expression matrix
       min_expression : float
           Minimum expression level
       min_samples : int
           Minimum samples with expression
       
       Returns
       -------
       pd.DataFrame
           Filtered expression matrix
       dict
           QC statistics
       """
       # Initial stats
       n_genes_initial = len(expr_matrix)
       n_samples_initial = len(expr_matrix.columns)
       
       # Filter low expression genes
       expressed = (expr_matrix > min_expression).sum(axis=1)
       keep_genes = expressed >= min_samples
       expr_filtered = expr_matrix[keep_genes]
       
       # Stats
       stats = {
           'initial_genes': n_genes_initial,
           'initial_samples': n_samples_initial,
           'filtered_genes': len(expr_filtered),
           'removed_genes': n_genes_initial - len(expr_filtered),
           'percent_kept': len(expr_filtered) / n_genes_initial * 100
       }
       
       return expr_filtered, stats

File I/O Utilities
------------------

Reading Expression Data
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def read_expression_matrix(file_path, log_transform=True, pseudocount=1):
       """
       Read and preprocess expression matrix.
       
       Parameters
       ----------
       file_path : str
           Path to expression file (CSV/TSV)
       log_transform : bool
           Apply log transformation
       pseudocount : float
           Pseudocount for log transformation
       
       Returns
       -------
       pd.DataFrame
           Expression matrix
       """
       # Detect separator
       with open(file_path, 'r') as f:
           first_line = f.readline()
           sep = '\t' if '\t' in first_line else ','
       
       # Read matrix
       expr = pd.read_csv(file_path, sep=sep, index_col=0)
       
       # Log transform if requested
       if log_transform:
           expr = np.log2(expr + pseudocount)
       
       return expr

Batch Processing
~~~~~~~~~~~~~~~~

.. code-block:: python

   def batch_process_species(input_dir, output_dir, process_func, **kwargs):
       """
       Process multiple species in batch.
       
       Parameters
       ----------
       input_dir : str
           Input directory with species subdirectories
       output_dir : str
           Output directory
       process_func : callable
           Function to apply to each species
       **kwargs
           Additional arguments for process_func
       
       Example
       -------
       def normalize_species(species_dir, output_dir):
           expr = pd.read_csv(species_dir / 'log_tpm.csv', index_col=0)
           expr_norm = normalize_expression(expr)
           expr_norm.to_csv(output_dir / 'log_tpm_norm.csv')
       
       batch_process_species('Input_Data', 'Output_Data', normalize_species)
       """
       from pathlib import Path
       import os
       
       input_path = Path(input_dir)
       output_path = Path(output_dir)
       output_path.mkdir(exist_ok=True)
       
       # Process each species
       for species_dir in input_path.iterdir():
           if not species_dir.is_dir():
               continue
           
           species_name = species_dir.name
           species_output = output_path / species_name
           species_output.mkdir(exist_ok=True)
           
           print(f"Processing {species_name}...")
           try:
               process_func(species_dir, species_output, **kwargs)
           except Exception as e:
               print(f"Error processing {species_name}: {e}")

Validation Utilities
--------------------

Data Consistency Checks
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def validate_multimodulon_data(mm):
       """
       Comprehensive validation of MultiModulon data.
       
       Parameters
       ----------
       mm : MultiModulon
           MultiModulon instance
       
       Returns
       -------
       dict
           Validation report
       """
       report = {
           'species': {},
           'alignment': {},
           'overall': {'valid': True, 'issues': []}
       }
       
       # Check each species
       for species in mm.species:
           species_report = {
               'has_expression': False,
               'has_metadata': False,
               'samples_match': False,
               'gene_table': False
           }
           
           data = mm[species]
           
           # Check expression data
           if data.log_tpm is not None:
               species_report['has_expression'] = True
               expr_samples = set(data.log_tpm.columns)
           else:
               report['overall']['issues'].append(f"{species}: No expression data")
               report['overall']['valid'] = False
           
           # Check metadata
           if data.sample_sheet is not None:
               species_report['has_metadata'] = True
               meta_samples = set(data.sample_sheet.index)
               
               # Check sample matching
               if species_report['has_expression']:
                   if expr_samples == meta_samples:
                       species_report['samples_match'] = True
                   else:
                       missing = expr_samples - meta_samples
                       extra = meta_samples - expr_samples
                       if missing:
                           report['overall']['issues'].append(
                               f"{species}: Samples in expression but not metadata: {missing}"
                           )
                       if extra:
                           report['overall']['issues'].append(
                               f"{species}: Samples in metadata but not expression: {extra}"
                           )
           
           # Check gene table
           if data.gene_table is not None:
               species_report['gene_table'] = True
           
           report['species'][species] = species_report
       
       # Check alignment
       if hasattr(mm, '_combined_gene_db') and mm._combined_gene_db is not None:
           report['alignment']['exists'] = True
           report['alignment']['species'] = list(mm._combined_gene_db.columns)
           report['alignment']['gene_families'] = len(mm._combined_gene_db)
       else:
           report['alignment']['exists'] = False
       
       return report

Usage Examples
--------------

Complete Utility Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from multimodulon import MultiModulon
   from multimodulon.gff_utils import gff2pandas, create_gene_table
   from multimodulon.utils.fasta_utils import extract_protein_sequences
   from pathlib import Path
   
   # Initialize
   mm = MultiModulon("Input_Data")
   
   # Step 1: Create gene tables from GFF
   create_gene_table(mm)
   
   # Step 2: Extract proteins if needed
   for species in mm.species:
       species_dir = Path(f"Input_Data/{species}")
       
       if not (species_dir / "protein.faa").exists():
           print(f"Extracting proteins for {species}")
           extract_protein_sequences(
               genome_fasta=species_dir / "genome.fna",
               gff_file=species_dir / "genome.gff",
               output_fasta=species_dir / "protein.faa"
           )
   
   # Step 3: Validate data
   validation_report = validate_multimodulon_data(mm)
   
   if not validation_report['overall']['valid']:
       print("Data issues found:")
       for issue in validation_report['overall']['issues']:
           print(f"  - {issue}")

Custom GFF Processing
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from multimodulon.gff_utils import gff2pandas
   
   # Parse multiple feature types
   gff_df = gff2pandas(
       'genome.gff',
       feature=['CDS', 'tRNA', 'rRNA', 'ncRNA']
   )
   
   # Analyze gene statistics
   feature_counts = gff_df['type'].value_counts()
   print("Feature statistics:")
   print(feature_counts)
   
   # Extract specific attributes
   if 'gene_biotype' in gff_df.columns:
       biotype_counts = gff_df['gene_biotype'].value_counts()
       print("\nGene biotypes:")
       print(biotype_counts)
   
   # Find overlapping genes
   for i, gene1 in gff_df.iterrows():
       overlaps = gff_df[
           (gff_df['seqid'] == gene1['seqid']) &
           (gff_df.index != i) &
           (gff_df['start'] < gene1['end']) &
           (gff_df['end'] > gene1['start'])
       ]
       if len(overlaps) > 0:
           print(f"{i} overlaps with: {overlaps.index.tolist()}")