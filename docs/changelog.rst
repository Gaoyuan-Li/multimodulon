Changelog
=========

Version 0.1.0 (2025-07-03)
-----------------------

Initial release of MultiModulon package.

**Features:**

* Multi-species/strain expression analysis using multi-view ICA
* Bidirectional Best Hits (BBH) for ortholog detection
* Union-Find based gene alignment across species
* Robust multi-view ICA with clustering
* Automated optimization of component numbers
* Comprehensive visualization functions
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

**Visualization:**

* Gene weight plots with COG coloring
* Activity plots with project/condition highlighting  
* Core component comparison across species
* Batch visualization support

**Optimization:**

* Single-gene filter–based optimization workflow
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
