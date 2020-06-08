#!/usr/bin/env python
"""Ensure that MuonGunGenerator DeepCore is working"""

import os
import tempfile
import shutil
from icecube import icetray, dataclasses, dataio, clsim
from I3Tray import I3Tray

try:

    from icecube.simprod.modules import MuonGunGenerator
    
    tmpdir = tempfile.mkdtemp(dir=os.getcwd())
    tmpfile = os.path.join(tmpdir,'test.i3')
    summaryfile = os.path.join(tmpdir,'summary.xml')
    gcdfile = os.path.expandvars('$I3_TESTDATA/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz')

    # prefer GPUs
    usegpus = any([device.gpu for device in clsim.I3CLSimOpenCLDevice.GetAllDevices()])    

    # make a very small DeepCore muon gun file
    n = MuonGunGenerator()
    n.SetParameter('outputfile',tmpfile)
    n.SetParameter('summaryfile',summaryfile)
    n.SetParameter('gcdfile',gcdfile)
    n.SetParameter('gamma',5.0)
    n.SetParameter('length',1600)
    n.SetParameter('radius',800)
    n.SetParameter('length_dc',500)
    n.SetParameter('radius_dc',150)
    n.SetParameter('deepcore',True)
    n.SetParameter('usegpus', usegpus)
    n.SetParameter('useonlydevicenumber',0)
    n.SetParameter('model','Hoerandel5_atmod12_SIBYLL')
    if n.Execute({}) != 0:
        raise Exception('MuonGunGenerator_DC did not return OK')

    # now check generated file
    tray = I3Tray()
    tray.Add('I3Reader', filename=tmpfile)
    def checky(frame):
        assert('I3MCTree' in frame)
        #print len(frame['I3MCTree'])
    tray.Add(checky, Streams=[icetray.I3Frame.DAQ])
    tray.Execute()
    
finally:
    shutil.rmtree(tmpdir)
