#!/usr/bin/env python

import sys
from I3Tray import I3Tray
from icecube import icetray
from icecube import dataio
from icecube import sim_services

if len(sys.argv) != 3: 
    print("usage: SplitMCTreeForHybridSimulation.py input_file output_file")
    sys.exit(-1)

tray = I3Tray()
tray.AddModule('I3Reader', Filename=sys.argv[1])
tray.AddModule('I3MCTreeHybridSimulationSplitter')
tray.AddModule('I3Writer', Filename=sys.argv[2])
tray.Execute()

