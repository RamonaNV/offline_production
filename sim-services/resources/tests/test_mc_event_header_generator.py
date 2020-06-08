#!/usr/bin/env python

import unittest

from I3Tray import I3Tray
from icecube import icetray
from icecube import dataclasses
from icecube import sim_services

class Source(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox("OutBox")

    def Configure(self):
        pass

    def Process(self):
        frame = icetray.I3Frame(icetray.I3Frame.DAQ)        
        self.PushFrame(frame)
        
class TestMCEventGenerator(unittest.TestCase):

    def setUp(self):
        self.tray = I3Tray()
        self.tray.AddModule(Source)

    def test_mc_event_generator_default(self):
                
        self.tray.AddModule("I3MCEventHeaderGenerator")                            

        def TestModule(frame):
            self.assertTrue("I3EventHeader" in frame)
            self.assertFalse(frame["MCTimeIncEventID"].value)
            
        self.tray.AddModule(TestModule, streams = [icetray.I3Frame.DAQ])

        self.tray.Execute(1)

    def test_mc_event_generator_increment_event_ID(self):
                
        self.tray.AddModule("I3MCEventHeaderGenerator",
                            IncrementEventID = True,
                            EventID = 100)
        self.expected_eid = 100
        def TestModule(frame):
            self.assertTrue("I3EventHeader" in frame)
            self.assertTrue(frame["MCTimeIncEventID"].value)
            self.assertEqual(frame["I3EventHeader"].event_id, self.expected_eid)
            self.expected_eid += 1
            
        self.tray.AddModule(TestModule, streams = [icetray.I3Frame.DAQ])

        self.tray.Execute(1000)

        
    def test_failure(self):
        # Can't specify both MJD and Year           
        self.tray.AddModule("I3MCEventHeaderGenerator",
                            MJD = 10333,
                            Year = 1999)
        self.assertRaises(RuntimeError, self.Execute)

    def Execute(self):
        self.tray.Execute(1)

if __name__ == '__main__':
    unittest.main()
