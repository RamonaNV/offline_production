.. _simclasses:

simclasses
==========

**Maintainer:** Alex Olivas - Though this is more of a community project.

This project is only meant to contain classes that inherit from I3FrameObject, or object that are stored in containers that inherit from I3FrameObject (e.g. I3Map or I3Vector).  The idea being that you should be able to include this in any other meta-project and be able to read data produced by simulation.  This keeps the dependencies minimal.

General
^^^^^^^
* :cpp:class:`I3MCPE`
* :cpp:class:`I3MCPulse`
* :cpp:class:`I3MMCTrack`

Photon Propagation
^^^^^^^^^^^^^^^^^^
* :cpp:class:`I3Photon`
* :cpp:class:`I3CompressedPhoton`
* :cpp:class:`I3ExtraGeometryItem`
* :cpp:class:`I3ExtraGeometryItemMove`
* :cpp:class:`I3ExtraGeometryItemUnion`
* :cpp:class:`I3ExtraGeometryItemCylinder`
* :cpp:class:`I3CylinderMap`

Shower Production
^^^^^^^^^^^^^^^^^
* :cpp:class:`I3CorsikaShowerInfo`
* :cpp:class:`I3MCNKGPoint`
* :cpp:class:`I3MCNKGInterpolation` 
* :cpp:class:`CorsikaLongStep`

BSM
^^^
* :cpp:class:`I3WimpParams` 


  Deprecated Classes
^^^^^^^^^^^^^^^^^^
The following classes are only kept around to read old data.
* :cpp:class:`I3GaussianPMTPulse`
* :cpp:class:`I3MCTWRParams`
* :cpp:class:`MMCWeight`
* :cpp:class:`I3MCPMTResponse`

Release Notes
^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 3
   
   release_notes
