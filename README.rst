.. image:: https://github.com/jdtibochab/coralme/blob/main/docs/logo.png

The **COmprehensive Reconstruction ALgorithm for ME-models (coralME)** is an automatic pipeline for the reconstruction of ME-models. coralME integrates existing ME-modeling packages `COBRAme`_, `ECOLIme`_, and `solveME`_, generalizes their functions for implementation on any prokaryote, and processes readily available organism-specific inputs for the automatic generation of a working ME-model.

coralME has four main objectives:

1. **Synchronize** input files to remove contradictory entries.
2. **Complement** input files from homology with a template organism to complete the E-matrix.
3. **Build** a working ME-model.
4. **Inform** the user about necessary steps to curate the ME-model.

For more information on coralME's inputs, arquitecture and usage, see coralME's `readthedocs`_.

Installation
------------

run ``pip install coralme``

Requirements
------------
coralME was tested with the following package versions:

- Python versions 3.7, 3.8, 3.9, or 3.10
- COBRApy version 0.25.0
- GUROBIpy version 9.5.2 (no license required)
- Ubuntu 22.04 is recommended (gfortran is required)
- Windows and MacOS users require to install `Gurobi`_ or IBM CPLEX Optimizer.

.. refs
.. _COBRAme: https://github.com/SBRG/cobrame
.. _ECOLIme: https://github.com/SBRG/ecolime
.. _solveME: https://github.com/SBRG/solvemepy
.. _readthedocs: https://coralme.readthedocs.io/
.. _Gurobi: https://www.gurobi.com/
.. _cplex: https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer
