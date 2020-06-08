#!/usr/bin/env python
import os
import subprocess
import threading
import time

from os.path import expandvars
from icecube import icetray, dataclasses,simclasses

import math



@icetray.traysegment
def MultiDomEffSample(tray, name,
    RandomService         = None,
    GeneratedEfficiency   = None,     # this should be the UnshadowedFraction used in clsim.
                                      #   If none, assmue the highest of the SampleEfficiencies
    SampleEfficiencies    = [0.9, 0.95, 0.99, 1.089,1.1979],  # Efficiencies to downsample to
    InputSeriesName       = "I3MCPESeriesMap",  # output series will be named [InputSeriesName]_[DomEff]
    DeleteOriginalSeries  = True,      # if true, remove original pulse series from frame after sample
    OverwriteOriginalSeries = False,   # useful if you only specify one SampleEfficiency
    ):
 
    """ This tray segment allows the downsampling of MCPE maps
    to perform mutliple dom efficiency simulations"""

    if GeneratedEfficiency is None:
        GeneratedEfficiency = max(SampleEfficiencies)

    for eff in SampleEfficiencies:
        outputname = InputSeriesName+"_"+str(eff)
        modname    = name+"_downsampleMCPE_"+str(eff)
        samplefraction = eff / GeneratedEfficiency
        tray.AddModule("I3DownsampleMCPE",modname,
                  RandomService = RandomService,
                  InputName     = InputSeriesName,
                  OutputName    = outputname,
                  SampleFrac    = samplefraction)

    if DeleteOriginalSeries:
        tray.AddModule('Delete',name+"_removeMCPEoriginals", Keys=[InputSeriesName])
    
    if OverwriteOriginalSeries and len(SampleEfficiencies) == 1:
        outputname = InputSeriesName+"_"+str(SampleEfficiencies[0])
        tray.AddModule('Rename',name+'_renameMCPE',
                       Keys=[outputname,InputSeriesName])
 
