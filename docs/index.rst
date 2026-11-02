
MBN-tools
=========

Python tools for running, processing, and analysing simulations
with MBN Explorer.

**MBN-tools** provides a collection of utilities for working with
MBN Explorer simulations, including simulation execution, file
handling, trajectory analysis, crystallography, and 3D visualisation.

Features
--------

* Run MBN Explorer simulation tasks directly from Python.
* Read, write, modify, and convert simulation files, including
  ``.task``, ``.xyz``, ``.dcd``, and ``.pdb``.
* Analyse trajectories using RDF, MSD, RMSD, and other methods.
* Perform crystallographic operations, including grouping atoms
  by crystallographic planes.
* Visualise atoms, planes, and simulation boxes in 3D using VisPy.

Installation
------------

Install MBN-tools using pip:

.. code-block:: bash

   pip install MBN-tools

Basic usage
-----------

.. code-block:: python

   import MBN_tools as MBN
   from MBN_tools import analysis, crystallography, visualise

   # Read an XYZ file
   xyz_data = MBN.read_xyz("simulation.xyz")

   # Calculate RMSD
   rmsd_values = analysis.calculate_rmsd(
       xyz_data,
       box_size=[50, 50, 100]
   )

   # Create a visualisation
   canvas = visualise.create_canvas()

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api

