.. _sim-services:

sim-services
============

**Maintainer** : Alex Olivas

.. toctree::
   :maxdepth: 3
   
   release_notes

Class and Namespace Overview
----------------------------

This project contains the following classes and namespaces.  

- :cpp:namespace:`I3SimConstants` - Namespace of various physics constants useful in simulation.
- :cpp:class:`I3PropagatorService` - Abstract base class for lepton propagator services.
- :cpp:class:`I3SumGenerator`  - A class for generating the sum of several random numbers which all have the same probability density.
- :cpp:class:`I3MCTreeHybridSimulationSplitter` - Splits an I3MCTree into two trees, one for tracks and one for cascades.
- :cpp:class:`I3MCPEtoI3MCHitConverter` - Converts an I3MCPESeriesMap to an I3MCHitSeriesMap.
- :cpp:class:`I3DownsampleMCPE` - Randomly downsample MCPEs from one collection to another.
- :cpp:class:`I3PropagatorModule` - Propagates all particles found in an MCTree.
- :cpp:class:`I3InIceCORSIKATrimmer` - Removes muons that have no chance of reaching the detector.
- :cpp:class:`I3CombineMCPE` - Combines several I3MCPEHitSeriesMaps into one.
- :cpp:class:`I3RemoveLargeDT` - Removes outlying I3MCPEs.

Deprecated Modules
------------------
* I3ModifyEventID - This modified more than the event ID *and* it was actually only used to set the start time, so this has been replaced with I3ModifyStartTime.
   
I3RemoveLargeDT 
---------------
.. toctree:: 
   :titlesonly: 

   remove_large_dt

I3PropagatorModule 
------------------
.. toctree:: 
   :titlesonly: 

   propagators

GCD Validation
----------------
.. toctree:: 
   :titlesonly: 

   gcd_validation

Code Review 
^^^^^^^^^^^ 
.. toctree:: 
   :titlesonly: 

   code_review
