#!/usr/bin/env python
"""Ensure that the NuGen API hasn't changed (too much)"""

import os
import tempfile
import shutil

from icecube.simprod.modules.nugen import NuGen 

from icecube import icetray, dataclasses, dataio
from I3Tray import I3Tray

try:
    tmpdir = tempfile.mkdtemp(dir=os.getcwd())
    tmpfile = os.path.join(tmpdir,'test.i3')
    summaryfile = os.path.join(tmpdir,'summary.xml')
    gcdfile = os.path.expandvars('$I3_TESTDATA/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz')
    
    # make a very small nugen file
    n = NuGen()
    n.SetParameter('nevents',1)
    n.SetParameter('outputfile',tmpfile)
    n.SetParameter('summaryfile',summaryfile)
    n.SetParameter('gcdfile',gcdfile)
    n.SetParameter('mjd',55697)
    n.SetParameter('NuFlavor','NuMu')
    if n.Execute({}) != 0:
        raise Exception('NuGen did not return OK')

    # now check generated file
    tray = I3Tray()
    tray.Add('I3Reader', filename=tmpfile)
    def checky(frame):
        assert('NuGPrimary' in frame)
        assert('I3MCTree_preMuonProp' in frame)
    tray.Add(checky, Streams=[icetray.I3Frame.DAQ])
    tray.Execute()

finally:
    shutil.rmtree(tmpdir)
