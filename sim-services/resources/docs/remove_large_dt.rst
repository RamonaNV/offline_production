Introduction
============
The goal of this module is to prevent isolated hits that are temporally separated from the rest
of the series from blowing up memory futher down the chain. Noise generator (Vuvuzela) will fill
the temporal space between with random noise. Similarly DOMLauncher will add beacon hits throughout
the time window. If the window is too large, this will blow up in memory. I3RemoveLargeDT calculates 
the median time of the I3MCPEs contained in the I3MCPESeries in each frame and will remove PEs that
are separated by more than half a time window MaxDeltaT which is configured as a module parameter.
The module thus eliminate hits that, while physically related to a triggered
event, would never be associated to that event by the DAQ.

Details
=======
The way I3RemoveLargeDT works is very simple. 

 #) Iterate through each I3MCPESeries in the I3MCPESeriesMap and sort if needed. (One can assume they are sorted already if they were genated by CLSim or PPC.)
 #) Build a vector of all PE times. At the same time. Keep track of the minimum and maximum time.
 #) If t_max - t_min < MaxDeltaT, no need to go any further since all hits are contained in the desired time window.
 #) Else: 

   #) Sort vector of PE times and calculate median time (middle element in vector).
   #) Iterate through I3MCPESeriesMap again and delete hits with time t, such that abs(t-t_med) > MaxDeltaT.



Usage
=====
The I3RemoveLargeDT module has three configurable parameters:
  * MaxDeltaT: Largest time span of PEs in an event. This is set by default to the size of the fixed-rate trigger.
  * InputResponse: Name of the input response series.
  * OutputResponse: Name of the output response series. Set different from InputResponse if you want to preserve original.
  * PreSorted: PEs are already sorted in time. This speeds up the process since it does not need to sort PEs to find earliest and latest hits.

This module should be added to the simulation chain after photon-propagation (e.g. PPC or Clsim) and before noise generation (Vuvuzela)
and PMTResponseSimulator.

.. code-block:: python

  import sim_services 

  tray.AddModule("I3RemoveLargeDT","cleanhits",
                 MaxDeltaT = 100*I3Units.ms,
                 InputResponse="I3MCPESeriesMap", 
                 OutputResponse="CleanedI3MCPESeriesMap", 
                 PreSorted=True)

