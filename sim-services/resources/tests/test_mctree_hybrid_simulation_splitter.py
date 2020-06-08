#!/usr/bin/env python

import sys
import unittest

from copy import deepcopy

from I3Tray import I3Tray
from icecube import icetray
from icecube import sim_services
from icecube import dataclasses
from icecube import phys_services

class Source(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox("OutBox")

    def Configure(self):
        pass

    def Process(self):
        frame = icetray.I3Frame(icetray.I3Frame.DAQ)        
        mctree = dataclasses.I3MCTree()

        cascade = dataclasses.I3Particle(dataclasses.I3Particle.Primary,
                                         dataclasses.I3Particle.EPlus)        
        muon = dataclasses.I3Particle(dataclasses.I3Particle.Primary,
                                      dataclasses.I3Particle.MuPlus)
        
        mctree.add_primary(cascade)
        mctree.add_primary(muon)
        
        frame["I3MCTree"] = mctree
        self.PushFrame(frame)

class TestI3MCTreeHybridSimulationSplitter(unittest.TestCase):

    def setUp(self):
        self.tray = I3Tray()
        self.tray.AddModule(Source)

    def test_down_sample_MCPE(self):
                        
        self.tray.AddModule("I3MCTreeHybridSimulationSplitter")

        def TestModule(frame):
            self.assertTrue("I3MCTreeTracks" in frame)
            self.assertTrue("I3MCTreeCascades" in frame)
            mctree_tracks = frame["I3MCTreeTracks"]
            mctree_cascades = frame["I3MCTreeCascades"]

            for particle in mctree_tracks :
                if particle.is_cascade:
                    self.assertEqual(particle.shape,
                                     dataclasses.I3Particle.Dark,
                                     "cascades in the track tree should be dark")
                else:
                    self.assertNotEqual(particle.shape,
                                     dataclasses.I3Particle.Dark,
                                     "only cascades should be dark in the track tree")

            for particle in mctree_cascades :
                if particle.is_track:
                    self.assertEqual(particle.shape,
                                     dataclasses.I3Particle.Dark,
                                     "tracks in the cascade tree should be dark")
                else:
                    self.assertNotEqual(particle.shape,
                                     dataclasses.I3Particle.Dark,
                                     "only tracks should be dark in the cascade tree")
                    
            
        self.tray.AddModule(TestModule, streams = [icetray.I3Frame.DAQ])
        self.tray.Execute(1)

    def Execute(self):
        self.tray.Execute(1)
        
if __name__ == '__main__':
    unittest.main()
