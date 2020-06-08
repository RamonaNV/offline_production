#!/usr/bin/env python
from I3Tray import *

from icecube import icetray, dataclasses, dataio, trigger_sim

import random

from os.path import expandvars
import sys

class TestModule(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('DOMsToExclude', 'These DOMs should be excluded.', [])

    def Configure(self):
        self.domsToExclude = self.GetParameter('DOMsToExclude')

    def DAQ(self, frame):
        print("TestModule : DAQ!!!")

            
        self.PushFrame(frame)
