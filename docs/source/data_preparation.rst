Data Preparation
================

This guide explains how to prepare your data for MultiModulon analysis.

Directory Structure
-------------------

MultiModulon expects data organized in a specific directory structure:

.. code-block:: text

   Input_Data/
   ├── Species1/
   │   ├── log_tpm.csv           # Expression matrix (required)
   │   ├── log_tpm_norm.csv      # Normalized expression (optional)
   │   ├── sample_sheet.csv      # Sample metadata (required)
   │   ├── gene_table.csv        # Gene annotations (optional)
   │   ├── genome.fna            # Genome sequence (for BBH)
   │   ├── genome.gff            # Gene annotations (for BBH)
   │   └── protein.faa           # Protein sequences (optional)
   ├── Species2/
   │   └── ... (same structure)
   └── Species3/
       └── ... (same structure)

Required Files
--------------

Expression Matrix (log_tpm.csv)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Format**: CSV file with genes as rows and samples as columns
* **Values**: Log-transformed TPM (Transcripts Per Million) values
* **Index**: Gene identifiers (must match gene_table if provided)

Example:

.. code-block:: text

   gene_id,Sample1,Sample2,Sample3
   gene001,5.2,4.8,5.1
   gene002,0.3,0.5,0.2
   gene003,7.1,7.3,6.9

Sample Sheet (sample_sheet.csv)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Format**: CSV file with samples as rows
* **Required columns**: None (index must match expression matrix columns)
* **Recommended columns**:
  
  - ``condition``: Experimental condition
  - ``project``: Project or study name
  - ``biological_replicate``: Replicate number (1, 2, 3, etc.)
  - ``study_accession``: Study identifier (e.g., from GEO)
  - ``sample_description``: Brief description

Example:

.. code-block:: text

   sample_id,condition,project,biological_replicate
   Sample1,Control,ProjectA,1
   Sample2,Control,ProjectA,2
   Sample3,Treatment,ProjectA,1

Optional Files
--------------

Gene Table (gene_table.csv)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Format**: CSV file with genes as rows
* **Index**: Must match expression matrix gene identifiers
* **Useful columns**:
  
  - ``gene_name``: Human-readable gene name
  - ``product``: Gene product description
  - ``COG``: COG category
  - ``start``, ``end``, ``strand``: Genomic coordinates

Normalized Expression (log_tpm_norm.csv)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Pre-normalized expression matrix
* Same format as log_tpm.csv
* If not provided, log_tpm will be used directly

Files for BBH Analysis
----------------------

To perform gene alignment across species, you need either:

Option 1: Genome Files
~~~~~~~~~~~~~~~~~~~~~~

* ``genome.fna``: Genome sequence in FASTA format
* ``genome.gff``: Gene annotations in GFF3 format

Option 2: Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``protein.faa``: Protein sequences in FASTA format

The protein sequences will be automatically extracted from genome files if not provided.

Data Quality Checklist
----------------------

Before running analysis, ensure:

1. ✓ All expression values are log-transformed
2. ✓ Sample names are consistent between expression matrix and sample sheet
3. ✓ Gene identifiers are consistent across all files
4. ✓ No missing values in expression matrix (or handle appropriately)
5. ✓ Biological replicates are properly labeled (1, 2, 3, not 1, 1, 2)
6. ✓ File encodings are UTF-8 (especially important for gene names)

Handling Special Cases
----------------------

Multiple Conditions
~~~~~~~~~~~~~~~~~~~

If you have multiple conditions in your experiment:

.. code-block:: python

   # The condition column will be automatically detected
   # for grouped visualization
   mm.view_iModulon_activities(
       species='Species1',
       component='Core_1',
       highlight_condition=['Treatment1', 'Treatment2']
   )

Missing Gene Annotations
~~~~~~~~~~~~~~~~~~~~~~~~

If gene annotations are not available:

.. code-block:: python

   # Create gene table from GFF file
   mm.create_gene_table()
   
   # Or use the gff2pandas utility
   gene_df = mm.gff2pandas('genome.gff', feature='CDS')

Large Datasets
~~~~~~~~~~~~~~

For datasets with many samples (>1000):

* Consider downsampling or splitting into batches
* Use GPU mode for faster computation
* Increase memory allocation if needed

Next Steps
----------

Once your data is properly formatted:

1. Follow the :doc:`quickstart` guide
2. See :doc:`initialization` for detailed initialization options
3. Check :doc:`gene_alignment` for BBH and alignment procedures