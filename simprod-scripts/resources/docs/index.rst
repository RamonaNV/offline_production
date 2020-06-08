.. highlight:: python

.. _simprod-scripts:

===============
Simprod-Scripts
===============

This project is a collection of scripts, tray segments and IceProd modules used in simulation 
production. The aim is to provide a central place with standard segments for running
simulation in both production and privately.

* IceProd modules: :doc:`ipmodule` are basic wrappers around tray segments that provide an interface for IceProd.
* Tray Segments: :doc:`segments` are IceTray meta-modules that contain several I3Modules with default parameters.
* Scripts: :doc:`scripts` is a collection of python scripts used in simulation production
* Examples: The directory simprod-scripts/resources/examples contains a collection of example scripts for running IPModules
* Tests

Rules and guidlines for scripts, segments, modules, and tests.

* Give one entity one cohesive responsibilty.
* Write documentation for the new IceCube member.  Make sure it's clear.  If you're having trouble writing clear documentation, perhaps you didn't give the entity you're writing documentation for one cohesive responsibility and maybe rethink your design.
* python/modules

  * One module per file - Don't include anything in here that IceProd does not call directly.

.. toctree::
    :maxdepth: 2
    :titlesonly:
    
    release_notes
    ipmodule
    modules
    segments
    scripts
    examples
    production_test
    corsika_tarballs
    
