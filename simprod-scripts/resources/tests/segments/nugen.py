#!/usr/bin/env python

# Ensure that the NuGen API hasn't changed

from icecube import icetray, dataclasses, dataio, phys_services
from icecube.simprod.segments.GenerateNeutrinos import GenerateNeutrinos
from I3Tray import I3Tray

tray = I3Tray()

tray.context['I3RandomService'] = phys_services.I3GSLRandomService(1337)
tray.Add("I3InfiniteSource")
tray.Add(GenerateNeutrinos, NumEvents=1)

def checky(frame):
    assert('NuGPrimary' in frame)
    assert('I3MCTree_preMuonProp' in frame)
tray.Add(checky, Streams=[icetray.I3Frame.DAQ])

tray.Execute()
