#!/usr/bin/env python

import unittest

from I3Tray import I3Tray
from icecube import icetray
from icecube import dataclasses
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
        mcpes = simclasses.I3MCPESeries()
        for i in range(1000) :            
            mcpes.append(simclasses.I3MCPE())
        mcpemap = simclasses.I3MCPESeriesMap()
        mcpemap[icetray.OMKey(21,30)] = mcpes
        frame["I3MCPESeriesMap"] = mcpemap

        self.PushFrame(frame)
        
class TestMCPEtoMCHitConverter(unittest.TestCase):

    def setUp(self):
        self.tray = I3Tray()
        self.tray.AddModule(Source)

    def test_mcpe_to_mchit_converter_default(self):
                
        self.tray.AddModule("I3MCPEtoI3MCHitConverter")

        def TestModule(frame):
            self.assertTrue("I3MCHitSeriesMap" in frame)
            self.assertTrue("I3MCPESeriesMap" in frame)
            self.assertEqual(len(frame["I3MCHitSeriesMap"]),
                             len(frame["I3MCPESeriesMap"]))

            n_mchits = sum([len(series) for omkey,series in frame["I3MCHitSeriesMap"].iteritems()])                           
            n_mcpes = sum([len(series) for omkey,series in frame["I3MCPESeriesMap"]])                           
            
            self.assertEqual(n_mchits, n_mcpes)                             

        self.tray.AddModule(TestModule, streams = [icetray.I3Frame.DAQ])

        self.tray.Execute(1)

    def test_failure(self):
        # no input found
        self.tray.AddModule("I3MCPEtoI3MCHitConverter", InputResponse = "foo")
        self.assertRaises(RuntimeError, self.Execute)

    def Execute(self):
        self.tray.Execute(1)

if __name__ == '__main__':
    unittest.main()
