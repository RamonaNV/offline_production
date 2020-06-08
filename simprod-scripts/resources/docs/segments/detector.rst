.. _DetectorSim:

DetectorSim
======================

Corresponding **IceProd module**: :ref:`detectors`

Detector simulation in IceCube
-------------------------------

The ``DetectorSim`` tray segment takes as input a GCD and MCPESeriesMap 
that may or may not cointain coincident MCPEs produced by `Polyplopia`_. 
If configured to do so, it will add noise MCPEs using `Vuvuzela`_.

It then calls the `DetectorResponse`_  segment from `DOMLauncher`_ which is responsible for
generating PMTResponse (`I3MCPulseSeriesMap`) and `DOMLaunchSeriesMap` objects. It also simulates the
triggers using the configuration from the GCD.

`Vuvuzela` takes input MCPEs and adds noise PEs resulting from thermal and correlated noise.

`PMTResponseSimulator` takes input MCPEs and:

- Adds a weight corresponding to the pulse charge that photon would yield.
- Generates pre-pulses, after-pulses, and late-pulses
- Applies time jitter
- Simulates saturation

`DOMLauncher` is responsible for simulating:

- Discriminator
- LC-logic
- Digitization
- Simulated effects
- Electronic noise in the digitizers
- Beacon launches (CPU triggered launches)
- The FPGA Clock phase
- RAPcal time uncertainty

`TriggerSim` is responsible for simulating all of the DAQ triggers.

API
---

.. autofunction:: icecube.simprod.segments.DetectorSim

