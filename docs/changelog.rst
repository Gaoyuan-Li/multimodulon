Changelog
=========

Version 1.0.0
-------------

First published version with full tutorial.

**Highlights:**

* End-to-end MultiModulon workflow tutorial
* Core component characterization tutorial
* Unique component characterization tutorial
* Multi-species/strain expression analysis with multi-view ICA
* BBH-based ortholog detection and gene alignment
* Core/unique component optimization with non-single-gene filtering
* Robust multi-view ICA and threshold optimization
* Visualization and interpretation utilities for iModulons

Version 0.2.0
-------------

**Visualization:**

* ``view_iModulon_weights`` - Gene weight plots per species
* ``view_iModulon_activities`` - Activity bar plots across samples
* ``view_iModulon_genes`` - Gene membership table per component
* ``view_core_iModulon_weights`` - Core component weights across species
* ``compare_core_iModulon`` - Gene membership comparison heatmap
* ``compare_core_iModulon_activity`` - Cross-species activity comparison
* ``plot_iM_conservation_bubble_matrix`` - Conservation bubble matrix summary
* ``show_iModulon_activity_change`` - Condition-vs-condition activity scatter
* ``show_gene_iModulon_correlation`` - Gene expression vs activity correlation
* ``core_iModulon_stability`` - Cross-species stability scores

Version 0.1.0
-------------

Initial release of MultiModulon package.

**Features:**

* Multi-species/strain expression analysis using multi-view ICA
* Bidirectional Best Hits (BBH) for ortholog detection
* Union-Find based gene alignment across species
* Robust multi-view ICA with clustering
* Automated optimization of component numbers
* GPU acceleration support
* JSON-based data persistence

**Core Functionality:**

* ``MultiModulon`` class for managing multi-species analysis
* ``SpeciesData`` container for individual species data
* BBH generation using containerized BLAST
* Gene alignment with customizable thresholds
* Multi-view ICA with PyTorch backend
* Single-gene component filtering
* Otsu's method for threshold optimization

**Optimization:**

* Single-gene filter-based optimization workflow
* Separate optimization for core and unique components
* Cross-validation support
* Automated parameter selection

**Utilities:**

* GFF parsing and gene table creation
* eggNOG annotation integration
* FASTA sequence extraction
* Data validation and quality control

**Documentation:**

* Comprehensive API documentation
* Step-by-step tutorials
* Example workflows
