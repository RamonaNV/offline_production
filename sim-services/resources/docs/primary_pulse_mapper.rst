Introduction
============
The goal of this module is to convert the I3ParticleIDMap which maps pulse indices to source particle IDs
and remap it to the root (primary) of the MCTree that the given particle belongs to. As a space-saving measure,
the propagated I3MCTree is often deleted and with it the particles that produced the I3MCPEs and corresponding I3MCPulses.
This propagated tree can be regenerated from the saved state of the RNG but the particle IDs will be different. 


Usage
=====
The I3PrimaryPulseMapper.cxx module has three configurable parameters:
  * I3MCTreeName: Name of I3MCTree to build map from 
  * InputParticleIDMapname: Name of the original ParticleI3Map
  * OutputParticleIDMapname: Name of the output ParticleI3Map 

This module should be added to the simulation chain at anypoint after PMTResponseSimulator as long 
as the original I3MCTree and ParticleI3Map have not been deleted.

.. code-block:: python

  import sim_services 

  tray.AddModule("I3PrimaryPulseMapper","mappulses")
                 I3MCTreeName="I3MCTree", 
                 InputParticleIDMapname="I3MCPulseSeriesMapParticleIDMap", 
                 OutputParticleIDMapname="I3MCPulseSeriesMapPrimaryIDMap") 

