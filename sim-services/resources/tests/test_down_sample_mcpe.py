#!/usr/bin/env python

import sys
import unittest

from copy import deepcopy

from I3Tray import I3Tray
from icecube import icetray
from icecube import simclasses
from icecube import sim_services
from icecube import phys_services

class Source(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox("OutBox")

    def Configure(self):
        pass

    def Process(self):
        frame = icetray.I3Frame(icetray.I3Frame.DAQ)        
        mcpes = simclasses.I3MCPESeries()
        for i in range(1000) :            
            mcpes.append(simclasses.I3MCPE(1))
        mcpemap = simclasses.I3MCPESeriesMap()
        mcpemap[icetray.OMKey(21,30)] = mcpes
        frame["I3MCPESeriesMap"] = mcpemap
        self.PushFrame(frame)

class TestI3DownsampleMCPE(unittest.TestCase):

    def setUp(self):
        self.tray = I3Tray()
        self.tray.AddModule(Source)

    def test_down_sample_MCPE(self):
                
        rng = phys_services.I3GSLRandomService(42)
        
        self.tray.AddModule("I3DownsampleMCPE",
                            SampleFrac = 0.50,
                            RandomService = rng)

        def TestModule(frame):
            self.assertTrue("DownsampledMCPEs" in frame)
            mcpemap = frame["DownsampledMCPEs"]
            npes = sum([1 for pe in mcpemap[icetray.OMKey(21,30)]])
            self.assertTrue(npes > 400)
            self.assertTrue(npes < 600)                    
            
        self.tray.AddModule(TestModule, streams = [icetray.I3Frame.DAQ])
        self.tray.Execute(1)

    def test_failure_no_rng(self):
        # user needs to specify the input.  there is no reasonable default.
        self.tray.AddModule("I3DownsampleMCPE")
        self.assertRaises(RuntimeError, self.Execute)

    def Execute(self):
        self.tray.Execute(1)
        
if __name__ == '__main__':
    unittest.main()
