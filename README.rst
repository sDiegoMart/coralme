.. image:: https://github.com/jdtibochab/coralme/blob/main/docs/logo.png

The **COmprehensive Reconstruction ALgorithm for ME-models (coralME)** is an automatic pipeline for the reconstruction of ME-models. coralME integrates existing ME-modeling packages `COBRAme`_, `ECOLIme`_, and `solveME`_, generalizes their functions for implementation on any prokaryote, and processes readily available organism-specific inputs for the automatic generation of a working ME-model.

coralME has four main objectives:

1. **Synchronize** input files to remove contradictory entries.
2. **Complement** input files from homology with a template organism to complete the E-matrix.
3. **Reconstruct** a ME-model.
4. **Troubleshoot** the ME-model to make it functional.

For more information on coralME's inputs, architecture and usage, see coralME's `readthedocs`_.

Installation
------------

run ``pip install coralme``

Requirements
------------

coralME was tested with the following package versions:

- Python3, versions 3.7, 3.8, 3.9, and 3.10
- COBRApy version 0.26.3
- GUROBIpy version 9.5.2 (license is required)
- Ubuntu 22.04 is recommended (libgfortran.so.5 is required to execute MINOS and quad MINOS)
- Windows and MacOS users require to install `Gurobi`_ or `IBM CPLEX Optimizer <cplex_>`_

Compiled MINOS and quad MINOS are provided here in ``*.so`` files under ``coralme/solver``, and have been compiled using:

- Python3, versions 3.7, 3.8, 3.9, and 3.10
- wheel 0.38.4
- numpy 1.21.6
- scipy 1.7.3
- cython 0.29.32
- cpython 0.0.6

The source code can be obtained from Prof. Michael A. Saunders at Stanford University

.. refs
.. _COBRAme: https://github.com/SBRG/cobrame
.. _ECOLIme: https://github.com/SBRG/ecolime
.. _solveME: https://github.com/SBRG/solvemepy
.. _readthedocs: https://coralme.readthedocs.io/
.. _Gurobi: https://www.gurobi.com/
.. _cplex: https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer
