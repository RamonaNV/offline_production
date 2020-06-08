#!/usr/bin/env python

import unittest

from copy import deepcopy

from I3Tray import I3Tray
from icecube import icetray
from icecube import simclasses
from icecube import sim_services

class Source(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox("OutBox")

    def Configure(self):
        pass

    def Process(self):
        frame = icetray.I3Frame(icetray.I3Frame.DAQ)        
        mcpe = simclasses.I3MCPE()
        mcpes = simclasses.I3MCPESeries()
        mcpes.append(mcpe)
        mcpemap0 = simclasses.I3MCPESeriesMap()
        mcpemap0[icetray.OMKey(21,30)] = mcpes
        mcpemap1 = simclasses.I3MCPESeriesMap()
        mcpemap1[icetray.OMKey(21,31)] = deepcopy(mcpes)
        frame["I3MCPESeriesMap0"] = mcpemap0
        frame["I3MCPESeriesMap1"] = mcpemap1
        self.PushFrame(frame)
        

class TestI3CombineMCPE(unittest.TestCase):

    def setUp(self):
        self.tray = I3Tray()
        self.tray.AddModule(Source)

    def test_combine_PE_maps(self):
                
        self.tray.AddModule("I3CombineMCPE",
                            InputResponses = ["I3MCPESeriesMap0","I3MCPESeriesMap1"])

        def TestModule(frame):
            self.assertTrue("CombinedMCPEs" in frame)
            self.assertEqual(len(frame["CombinedMCPEs"]), 2)
            
        self.tray.AddModule(TestModule, streams = [icetray.I3Frame.DAQ])
        self.tray.Execute(1)

    def test_combine_PE_maps_user_defined_output(self):
                
        self.tray.AddModule("I3CombineMCPE",
                            InputResponses = ["I3MCPESeriesMap0","I3MCPESeriesMap1"],
                            OutputResponse = "NewOutput")

        def TestModule(frame):
            self.assertTrue("NewOutput" in frame)
            self.assertEqual(len(frame["NewOutput"]), 2)
            
        self.tray.AddModule(TestModule, streams = [icetray.I3Frame.DAQ])
        self.tray.Execute(1)

    def test_failure(self):
        # user needs to specify the input.  there is no reasonable default.
        self.tray.AddModule("I3CombineMCPE")
        self.assertRaises(RuntimeError, self.Execute)

    def Execute(self):
        self.tray.Execute(1)
        
if __name__ == '__main__':
    unittest.main()
