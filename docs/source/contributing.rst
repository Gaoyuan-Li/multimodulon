Contributing
============

We welcome contributions to the MultiModulon project! This guide will help you get started.

How to Contribute
-----------------

Reporting Issues
~~~~~~~~~~~~~~~~

If you find a bug or have a feature request:

1. Check if the issue already exists
2. Create a new issue with a clear title
3. Provide a minimal reproducible example
4. Include your environment details (OS, Python version, package versions)

Submitting Pull Requests
~~~~~~~~~~~~~~~~~~~~~~~~

1. Fork the repository
2. Create a feature branch (``git checkout -b feature-name``)
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request

Development Setup
-----------------

Setting Up Environment
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/yourusername/multimodulon.git
   cd multimodulon
   
   # Create virtual environment
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   
   # Install in development mode
   pip install -e .
   pip install -r requirements-dev.txt

Running Tests
~~~~~~~~~~~~~

.. code-block:: bash

   # Run all tests
   pytest
   
   # Run with coverage
   pytest --cov=multimodulon
   
   # Run specific test file
   pytest tests/test_core.py

Code Style
----------

We follow PEP 8 guidelines with these additions:

* Maximum line length: 100 characters
* Use type hints where appropriate
* Add docstrings to all public functions

.. code-block:: python

   def calculate_something(data: pd.DataFrame, threshold: float = 0.5) -> float:
       """
       Calculate something from the data.
       
       Parameters
       ----------
       data : pd.DataFrame
           Input data
       threshold : float, optional
           Threshold value (default: 0.5)
           
       Returns
       -------
       float
           Calculated result
       """
       # Implementation here
       return result

Documentation
-------------

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cd docs
   make html
   # View at docs/build/html/index.html

Writing Documentation
~~~~~~~~~~~~~~~~~~~~~

* Use NumPy style docstrings
* Include examples in docstrings
* Update relevant .rst files for new features
* Add new functions to API reference

Testing Guidelines
------------------

Test Structure
~~~~~~~~~~~~~~

.. code-block:: python

   import pytest
   from multimodulon import MultiModulon
   
   class TestMultiModulon:
       @pytest.fixture
       def sample_data(self):
           """Provide sample data for tests."""
           return create_sample_data()
       
       def test_initialization(self, sample_data):
           """Test MultiModulon initialization."""
           mm = MultiModulon(sample_data)
           assert len(mm.species) > 0
       
       def test_invalid_input(self):
           """Test handling of invalid input."""
           with pytest.raises(ValueError):
               MultiModulon("nonexistent_path")

Test Coverage
~~~~~~~~~~~~~

* Aim for >80% code coverage
* Test edge cases and error conditions
* Include integration tests for workflows
* Test both CPU and GPU modes (if applicable)

Areas for Contribution
----------------------

Current Priorities
~~~~~~~~~~~~~~~~~~

* **Performance optimization**: Improve speed for large datasets
* **New algorithms**: Additional ICA variants, optimization metrics
* **Visualization**: Interactive plots, new plot types
* **Documentation**: Tutorials, use cases, scientific background
* **Testing**: Increase coverage, add benchmarks

Feature Ideas
~~~~~~~~~~~~~

* Time-series analysis support
* Integration with other bioinformatics tools
* Web interface for visualization
* Automated biological interpretation
* Cross-species network inference

Review Process
--------------

Pull Request Review
~~~~~~~~~~~~~~~~~~~

1. Automated tests must pass
2. Code review by maintainer
3. Documentation updated if needed
4. Changelog entry added
5. Merge after approval

Review Criteria
~~~~~~~~~~~~~~~

* Code quality and style
* Test coverage
* Documentation completeness
* Performance impact
* Backward compatibility

Community
---------

Getting Help
~~~~~~~~~~~~

* GitHub Issues for bugs/features
* Discussions for questions
* Email maintainers for sensitive issues

Code of Conduct
~~~~~~~~~~~~~~~

* Be respectful and inclusive
* Welcome newcomers
* Focus on constructive feedback
* Collaborate openly

Recognition
~~~~~~~~~~~

Contributors will be:

* Listed in AUTHORS file
* Acknowledged in changelog
* Credited in relevant documentation

Thank you for contributing to MultiModulon!