Changelog
=========

Version 0.1.0 (2024-01)
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
* Cohen's d effect size filtering
* Otsu's method for threshold optimization

**Visualization:**

* Gene weight plots with COG coloring
* Activity plots with project/condition highlighting  
* Core component comparison across species
* Batch visualization support

**Optimization:**

* NRE and effect size based optimization
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
* Best practices guide

Known Issues
------------

* Font warnings on some systems - specify font_path parameter
* Memory usage high for >10 species - use batch processing
* GPU memory limits for very large datasets - reduce batch size

Future Plans
------------

* Support for time-series analysis
* Integration with regulatory databases
* Interactive visualization dashboard
* Automated biological interpretation
* Cross-species regulatory network inference