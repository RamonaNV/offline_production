import os
import sys
import getopt
import logging
from os.path import expandvars

from I3Tray import *
from icecube import icetray, dataclasses, simclasses, dataio
from icecube.icetray import I3Frame
from .. import ipmodule
from .. import segments

logger = logging.getLogger('polyplopia')



class PolyplopiaModule(ipmodule.ParsingModule):
    """
    Add background coincidences to signal MC
    """

    def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('gcdfile', 'Geometry Calibration Det-Status file', '')
        self.AddParameter('inputfile', 'Input filename', '')
        self.AddParameter('backgroundfile', 'Background filename', '')
        self.AddParameter('backgroundrate', 'Background rate (Hz) (don\'t use with corsika)', float('nan'))
        self.AddParameter('mctype', 'Type of primary simulation', 'corsika')
        self.AddParameter('outputfile', 'Output filename', 'output.i3.gz')
        self.AddParameter('seed', 'RNG Seed', 123)
        self.AddParameter('procnum', 'RNG stream number', 0)
        self.AddParameter('nproc', 'RNG number of streams', 1)
        self.AddParameter('PROPOSALParams','any other parameters for proposal',dict())
        self.AddParameter('Propagate','propagete particles and photons',True)
        self.AddParameter("GPU", 
                          "Graphics Processing Unit number (shoud default to environment if None)",
                          None)
        self.AddParameter("UseGPUs", 
                          "Use Graphics Processing Unit for photon propagation.", True)
        self.AddParameter("UsePPC", 
                          "Use PPC for photon propagation instead of CLSim.", False)
        self.AddParameter("oversize","over-R: DOM radius oversize scaling factor",5)
        self.AddParameter("holeiceparametrization", "Location of hole ice param files",
                                                  expandvars("$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm"))
        self.AddParameter("totalenergytoprocess", "Accumulate frames according to amount of light deposited in the detector", 0)
        self.AddParameter("efficiency","overall DOM efficiency correction",[0.99])
        self.AddParameter("volumecyl","set volume to regular cylinder (set to False for 300m spacing from the DOMs)",True)
        self.AddParameter("IceModelLocation","Location of ice model param files", expandvars("$I3_BUILD/ice-models/resources/models")) 
        self.AddParameter("IceModel","ice model subdirectory", "spice_3.2") 
        self.AddParameter("MCTreeName","Name of MCTree frame object", "I3MCTree") 
        self.AddParameter("PhotonSeriesName","Photon Series Name","I3MCPESeriesMap") 
        self.AddParameter("RawPhotonSeriesName","Raw Photon Series Name",None) 
        self.AddParameter("UseGeant4","Enable Geant4 propagation",False) 
        self.AddParameter('HistogramFilename', 'Histogram filename.', None)
        self.AddParameter('TimeWindow', 'Coincidence time window', 40.*I3Units.microsecond)
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 
        
    def Execute(self, stats):
        'Runs the tray, but the important stuff happens in segment_X(...)'

        if not ipmodule.ParsingModule.Execute(self,stats):
            return 0

        # Load libraries
        from icecube import (earthmodel_service, PROPOSAL, cmc, phys_services)
        from icecube import (polyplopia, dataclasses, dataio)

        # support json ParamsMap, so we can get our dict from the iceprod config file
        try:
            import json
        except ImportError:
            json = None

        # Instantiate a tray
        tray = I3Tray()
        randomService = phys_services.I3SPRNGRandomService(
             seed = self.seed, 
             nstreams = self.nproc, 
             streamnum = self.procnum) \
            if not self.usegslrng else phys_services.I3GSLRandomService(self.seed)

        tray.context['I3RandomService'] = randomService

        # Configure IceTray modules
        tray.AddModule("I3Reader","reader", FilenameList=[self.gcdfile,self.inputfile])

        from .. import segments
        from I3Tray import I3Units

        tray.AddSegment(segments.PolyplopiaPhotons,"coincifypes",
               RandomService=randomService,
               mctype=self.mctype,
               mctree_name = "I3MCTree",
               bgfile = self.backgroundfile,
               GCDFile = self.gcdfile,
               timewindow = self.timewindow,
               GPU = self.gpu,
               UseGPUs = self.usegpus,
               UsePPC = self.useppc,
               IceModel = self.icemodel,
               IceModelLocation = self.icemodellocation,
               DOMOversizeFactor = self.oversize,
               HoleIceParameterization = self.holeiceparametrization,
               Efficiency = self.efficiency,
               PhotonSeriesName = self.photonseriesname,
               PROPOSALParams=self.proposalparams,
               rate = self.backgroundrate*I3Units.hertz
           ) 

        if self.histogramfilename:         
            from icecube.production_histograms import ProductionHistogramModule
            from icecube.production_histograms.histogram_modules.simulation.mcpe_module import I3MCPEModule
        
            tray.AddModule(ProductionHistogramModule, 
                           Histograms = [I3MCPEModule],
                           OutputFilename = self.histogramfilename)

        tray.AddModule("I3Writer","writer",
            Filename=self.outputfile,
            Streams=[icetray.I3Frame.TrayInfo,
                     icetray.I3Frame.DAQ,
                     icetray.I3Frame.Stream('S'),
                     icetray.I3Frame.Stream('M')])
         

        # Execute the Tray
        tray.Execute()        

        # Free memory
        del tray
        return 0


class PolyplopiaMCPEMerge(ipmodule.ParsingModule):
    """
    Add background coincidences to signal MC
    """

    def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('gcdfile', 'Geometry Calibration Det-Status file', '')
        self.AddParameter('inputfile', 'Input filename', '')
        self.AddParameter('backgroundfile', 'Background filename', '')
        self.AddParameter('mctype', 'Type of primary simulation', 'corsika')
        self.AddParameter('outputfile', 'Output filename', 'output.i3.gz')
        self.AddParameter('seed', 'RNG Seed', 123)
        self.AddParameter('procnum', 'RNG stream number', 0)
        self.AddParameter('nproc', 'RNG number of streams', 1)
        self.AddParameter('TimeWindow', 'Coincidence time window', 40.*I3Units.microsecond)
        self.AddParameter("MCTreeName","Name of MCTree frame object", "I3MCTree") 
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 

    def Execute(self, stats):
        'Runs the tray, but the important stuff happens in segment_X(...)'

        if not ipmodule.ParsingModule.Execute(self,stats):
            return 0

        # Load libraries
        from icecube import (earthmodel_service, PROPOSAL, cmc, phys_services)
        from icecube import (polyplopia, dataclasses, dataio)

        # support json ParamsMap, so we can get our dict from the iceprod config file
        try:
            import json
        except ImportError:
            json = None

        # Instantiate a tray
        tray = I3Tray()

        randomService = phys_services.I3SPRNGRandomService(
             seed = self.seed, 
             nstreams = self.nproc, 
             streamnum = self.procnum)\
            if not self.usegslrng else phys_services.I3GSLRandomService(seed)

        tray.context['I3RandomService'] = randomService

        # Configure IceTray modules
        tray.AddModule("I3Reader","reader", FilenameList=[self.gcdfile,self.inputfile])

        from .. import segments
        from I3Tray import I3Units

        tray.Add(segments.PolyplopiaMergePEs,"merge",
               mctype=self.mctype,
               RandomService=randomService,
               mctree_name = "I3MCTree",
               bgfile = self.backgroundfile,
               separate_coincident_mctree_name = "BackgroundI3MCTree",
               timewindow = self.timewindow,
               rate = 5.0*I3Units.kilohertz)


        tray.AddModule("I3Writer","writer",
            Filename=self.outputfile,
            Streams=[icetray.I3Frame.TrayInfo,
                     icetray.I3Frame.DAQ,
                     icetray.I3Frame.Stream('S'),
                     icetray.I3Frame.Stream('M')])

        # Execute the Tray
        tray.Execute()
        

        # Free memory
        del tray
        return 0

