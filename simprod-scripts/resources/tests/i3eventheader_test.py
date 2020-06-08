#!/usr/bin/env python
"""
Verify the length of the I3EventHeader, compared to triggers/hits/pulses.
"""

import os

from I3Tray import I3Tray
from icecube import icetray, dataio, dataclasses, phys_services

from icecube.simprod.segments import GenerateNeutrinos, PropagateMuons, PropagatePhotons, DetectorSim

gcdfile = os.path.expandvars(os.path.expanduser('$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'))

randomService = phys_services.I3GSLRandomService(12345)

tray = I3Tray()
tray.context['I3RandomService'] = randomService
tray.Add("I3InfiniteSource",
         Prefix=gcdfile,
         Stream=icetray.I3Frame.DAQ
)
tray.Add(GenerateNeutrinos,
         NumEvents=10,
         FromEnergy=10000,
         ToEnergy=50000,
         RandomService=randomService,
)
tray.Add(PropagateMuons,
         RandomService=randomService,
)
tray.Add(PropagatePhotons, UseAllCPUCores=True)
tray.Add(DetectorSim,
         RunID=1,
         GCDFile=gcdfile,
         RandomService='I3RandomService',
)

def verifyheader(fr):
    h = fr['I3EventHeader']
    header_time = h.end_time-h.start_time+0.1 # extra 0.1 because truncation
    trigger_time = max(t.time+t.length for t in fr['I3TriggerHierarchy'] if t.fired)
    if trigger_time > header_time:
        print('Bad Header')
        print('header duration:',header_time)
        print('trigger duration:',trigger_time)
        print(h)
        raise Exception('Bad Header')
tray.Add(verifyheader, Streams=[icetray.I3Frame.DAQ])
tray.Execute()
