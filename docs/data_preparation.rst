Data Preparation
================

This guide explains how to prepare your data for MultiModulon analysis.

Directory Structure
-------------------

MultiModulon expects data organized in a specific directory structure (Output from https://github.com/Gaoyuan-Li/MAPPED):

.. code-block:: text

   Input_Data/
   ├── Species1/
   │   ├── samplesheet/
   │   │   ├── sample_sheet.csv      # Sample metadata (required)
   │   ├── expression_matrices/
   │   │   ├── log_tpm.csv           # Expression matrix (required)
   │   │   └── log_tpm_norm.csv      # Normalized expression (required)
   │   │   └── counts.csv            # Counts matrix (optional)
   │   │   └── tpm.csv               # TPM matrix (optional)
   │   ├── ref_genome/
   │   │   ├── genome.fna            # Genome sequence (required)
   │   │   ├── genome.gff            # Gene annotations (required)
   │   │   └── protein.faa           # Protein sequences (required)
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
  
  - ``condition``: Experimental condition # only when available
  - ``project``: Project or study name # only when available
  - ``biological_replicate``: Replicate number (1, 2, 3, etc.) # only when available
  - ``study_accession``: Study identifier (e.g., from GEO) # only when available
  - ``sample_description``: Brief description # only when available

Example:

.. code-block:: text

   sample_id,condition,project,biological_replicate
   Sample1,Control,ProjectA,1
   Sample2,Control,ProjectA,2
   Sample3,Treatment,ProjectA,1

Normalized Expression (log_tpm_norm.csv)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Pre-normalized expression matrix
* Same format as log_tpm.csv
* If not provided, log_tpm will be used directly

Files for BBH Analysis and Gene Annotation
----------------------

To perform gene alignment across species and add gene annotation, you need:

* ``genome.fna``: Genome sequence in FASTA format
* ``genome.gff``: Gene annotations in GFF3 format
* ``protein.faa``: Protein sequences in FASTA format