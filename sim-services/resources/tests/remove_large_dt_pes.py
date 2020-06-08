#!/usr/bin/env python
import random
import sys

from icecube import icetray
from icecube import simclasses
from icecube import dataio
from icecube import sim_services

from I3Tray import I3Tray
from I3Tray import I3Units

class pe_inserter(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("TimeWindow", "Tuple of min and max times", None)
        self.AddOutBox("OutBox")
        
    def Configure(self) :
        self.time_window = self.GetParameter("TimeWindow")
        
    def DAQ(self, frame):
        min_time = self.time_window[0]
        max_time = self.time_window[1]

        peseries = simclasses.I3MCPESeries()
        
        # put 1000 within the time window
        for i in range(1000):
            pe = simclasses.I3MCPE()
            pe.time = random.uniform(min_time/10., max_time/10.)
            peseries.append(pe)

        # determine the median time
        times = sorted([pe.time for pe in peseries])
        median_time = times[int(len(times)/2.)]
                   
        # put 250 outside the time window
        for i in range(250):
            pe = simclasses.I3MCPE()
            time = min_time - random.uniform(100*I3Units.ms, 1000*I3Units.ms)
            pe.time = time
            peseries.append(pe)

        for i in range(250):
            pe = simclasses.I3MCPE()
            time = max_time + random.uniform(100*I3Units.ms, 1000*I3Units.ms)
            pe.time = time
            peseries.append(pe)

        pemap = simclasses.I3MCPESeriesMap()
        pemap[icetray.OMKey(21,30)] = peseries

        frame["I3MCPESeriesMap"] = pemap
        
        self.PushFrame(frame)

def check_cleaning(frame):
    pemap = frame["CleanedI3MCPESeriesMap"]

    npe = len([1 for omkey, peseries in pemap for pe in peseries])

    if npe != 1000 :
        print("There should be 1000 PEs.  Got %d." % npe)
        sys.exit(1)


tray = I3Tray()

tray.Add("I3InfiniteSource")

tray.AddModule(pe_inserter, \
               TimeWindow = (-50*I3Units.ms, 50*I3Units.ms))

tray.AddModule("I3RemoveLargeDT", PreSorted = False)

tray.AddModule(check_cleaning, Streams = [icetray.I3Frame.DAQ])
               
tray.Execute(1)

