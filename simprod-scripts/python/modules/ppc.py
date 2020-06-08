#!/usr/bin/env python
################################ Tray 1 ##############################
#   
#   IceProd classes for processing hit-generation in IceSim.
#   Author: juancarlos@icecube.wisc.edu
#   
######################################################################
import os
import tempfile

from os.path import expandvars, join

from I3Tray import *
from .. import ipmodule
from .clsim import save_stats
from ..util import ReadI3Summary, WriteI3Summary


def PPCTraySegment(tray,
                   name,
                   UnshadowedFraction = 1.,
                   DOMOversizeFactor = 1,
                   gpulib = "opencl",
                   volumecyl = True,
                   IceModelLocation = expandvars("$I3_SRC/ice-models/resources/models"),
                   keep_empty_events=False,
                   HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm"),
                   IceModel="Spice3.2.1",
                   InputMCTree="I3MCTree",
                   MCPESeriesName = "I3MCPESeriesMap",
                   UseGPUs=True,
                   GPU=-1,
                   tempdir = None):

        """
        PPC Photon Propagation Code TraySegment (supports CUDA/OpenCL) 
        """

        # Load libraries 
        from icecube import icetray, interfaces ,dataclasses ,simclasses,sim_services
        from icecube import polyplopia
        from icecube import ppc
        
        # Do one or the other
        UseCPU = not UseGPUs

        ppcIceModel = None
        if IceModelLocation is None:
             IceModelLocation = expandvars("$I3_SRC/ice-models/resources/models")

        if IceModel == "SpiceMie":
             ppcIceModel = expandvars(IceModelLocation+"/spice_mie")
        elif IceModel == "SpiceLea":
             ppcIceModel = expandvars(IceModelLocation+"/spice_lea")
        elif IceModel == "Spice3":
             ppcIceModel = expandvars(IceModelLocation+"/spice_3")
        elif IceModel == "Spice3.1":
             ppcIceModel = expandvars(IceModelLocation+"/spice_3.1")
        elif IceModel == "Spice3.2":
             ppcIceModel = expandvars(IceModelLocation+"/spice_3.2")
        elif IceModel == "Spice3.2.1":
             ppcIceModel = expandvars(IceModelLocation+"/spice_3.2.1")
        elif os.path.exists(expandvars(IceModelLocation+"/"+IceModel)):
             ppcIceModel = expandvars(IceModelLocation+"/"+IceModel)
        else:
             raise RuntimeError("Unknown ice model: %s", IceModel)
                    
        os.putenv("ICEMODELDIR", ppcIceModel)
        os.putenv("PPCHOLEICE", HoleIceParameterization)

        if UseGPUs:
            os.putenv("OGPU","1")
        else:
            os.putenv("OCPU","1")
        if GPU >= 0 and UseGPUs:
            os.putenv("CUDA_VISIBLE_DEVICES",str(GPU))
            os.putenv("COMPUTE",":0."+str(GPU))
            os.putenv("GPU_DEVICE_ORDINAL",str(GPU))

        tray.AddModule("i3ppc", "ppc", 
                  If = lambda f: f[InputMCTree].size() or keep_empty_events,
                  gpu=GPU, 
                  efficiency_scaling_factor = UnshadowedFraction,
                  cyl=volumecyl,
                  keep=keep_empty_events,
                  MCTree=InputMCTree)

        # PPC does not have an option for setting the name of the PE map.
        # If the default name of PE map changes, this needs to be updated.
        from icecube.simclasses import I3MCPESeriesMap
        def add_empty_pes(f,MCPESeriesName="MCPESeriesMap"):
                  if MCPESeriesName not in f:
                     f[MCPESeriesName] = I3MCPESeriesMap()
        tray.Add(add_empty_pes,MCPESeriesName="MCPESeriesMap",streams=[icetray.I3Frame.DAQ])
        if MCPESeriesName != "MCPESeriesMap":
           tray.AddModule("Rename", keys=["MCPESeriesMap",MCPESeriesName])


###### IP Modules ###########

class PPC(ipmodule.ParsingModule):
   """
    GPU Photon propagation
   """
   def __init__(self):
        ipmodule.ParsingModule.__init__(self)
        self.configured = False

        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('inputfilelist','list of input filenames',[])
        self.AddParameter('outputfile','Output filename','')
        self.AddParameter('seed','RNG seed',0)
        self.AddParameter('procnum','job number (RNG stream number)',0)
        self.AddParameter('nproc','Number of jobs (Number of RNG streams)',1)
        self.AddParameter('summaryfile','JSON Summary filename','summary.json')
        self.AddParameter('MJD','MJD (0=do not modify)',0)
        self.AddParameter("GPU", 
                          "Graphics Processing Unit number (shoud default to environment if None)",
                          -1)
        self.AddParameter("UseGPUs", "Use Graphics Processing Unit",True)
        self.AddParameter('RunMPHitFilter',"Run polyplopia's mphitfilter",True)
        self.AddParameter("oversize","over-R: DOM radius oversize scaling factor",5)
        self.AddParameter("efficiency","overall DOM efficiency scaling factor (systematics)",1.00)
        self.AddParameter("gpulib","set gpu library to load (defaults to cuda)","opencl")
        self.AddParameter("volumecyl","set volume to regular cylinder (set to False for 300m spacing from the DOMs)",True)
        self.AddParameter("PhotonSeriesName","Photon Series Name","I3MCPESeriesMap") 
        self.AddParameter("IceModelLocation","Location of ice model param files", expandvars("$I3_BUILD/ice-models/resources/models")) 
        self.AddParameter("IceModel","ice model subdirectory", "SpiceMie") 
        self.AddParameter("holeiceparametrization", "Location of hole ice param files", 
                          expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.h2-50cm"))
        self.AddParameter("MCTreeName","Name of MCTree frame object", "I3MCTree") 
        self.AddParameter("KeepEmptyEvents","Don't discard events with no MCPEs", False) 
        self.AddParameter('HistogramFilename', 'Histogram filename.', None)
        self.AddParameter('EnableHistogram', 'Write a SanityChecker histogram file.', False)
        self.AddParameter('PropagateMuons', 'Run PROPOSAL to do in-ice propagation', True)
        self.AddParameter('PROPOSALParams','any other parameters for proposal',dict())
        self.AddParameter('TempDir', 'Temporary working directory with the ice model', None)
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 
        self.configured = False

   def Configure(self,tray):
        from icecube import icetray, dataclasses, dataio, phys_services, interfaces
        from .ppc import PPCTraySegment
        from I3Tray import I3Tray,I3Units


        from .. import segments
        if self.propagatemuons:
        	randomServiceForPropagators = phys_services.I3SPRNGRandomService(
             		seed = self.seed,
             		nstreams = self.nproc*2,
             		streamnum = self.nproc + self.procnum)\
                if not self.usegslrng else phys_services.I3GSLRandomService(seed = self.seed*self.nproc+self.procnum)

        	tray.context['I3PropagatorRandomService'] = randomServiceForPropagators

        	tray.AddModule("Rename","rename_corsika_mctree",Keys=['I3MCTree','I3MCTree_preMuonProp'])
        	tray.AddSegment(segments.PropagateMuons, 'propagator',
                        RandomService= randomServiceForPropagators,
                        **self.proposalparams
        	) 


        
        tray.AddSegment(PPCTraySegment,"ppc_photons",
			gpu = self.gpu,
			usegpus = self.usegpus,
			UnshadowedFraction = self.efficiency,
			DOMOversizeFactor = self.oversize,
			IceModelLocation = self.icemodellocation,
			HoleIceParameterization = self.holeiceparametrization,
			IceModel = self.icemodel,
			volumecyl = self.volumecyl,
			gpulib = self.gpulib,
			InputMCTree = self.mctreename,
			keep_empty_events = self.keepemptyevents,
			MCPESeriesName = self.photonseriesname,
                        tempdir = self.tempdir)

        if self.runmphitfilter:
           tray.AddModule("MPHitFilter","hitfilter",
              HitOMThreshold=1,
              RemoveBackgroundOnly=False,
              I3MCPESeriesMapName=self.photonseriesname)

        if self.enablehistogram and self.histogramfilename:         
            from icecube.production_histograms import ProductionHistogramModule
            from icecube.production_histograms.histogram_modules.simulation.mcpe_module import I3MCPEModule
        
            tray.AddModule(ProductionHistogramModule, 
                           Histograms = [I3MCPEModule],
                           OutputFilename = self.histogramfilename)


        tray.AddModule("I3Writer","writer", 
            filename=self.outputfile,
            Streams=[icetray.I3Frame.TrayInfo,
                     icetray.I3Frame.DAQ,
                     icetray.I3Frame.Stream('S'),
                     icetray.I3Frame.Stream('M')])

   def Execute(self,stats):

        from icecube import dataclasses

        if not ipmodule.ParsingModule.Execute(self,stats):
                return 0

        # Instantiate a tray 
        tray = I3Tray()

        # Configure IceTray services
        if os.path.exists(self.summaryfile):
           summary = ReadI3Summary(self.summaryfile)
        else: 
           summary = dataclasses.I3MapStringDouble()
        tray.context['I3SummaryService'] = summary

        # Configure IceTray modules 
        tray.AddModule("I3Reader", "reader",filenamelist=[self.gcdfile]+self.inputfilelist)

        self.Configure(tray)

        # Execute the Tray
        tray.Execute()

        summary = tray.context['I3SummaryService']
        save_stats(tray,summary,stats)
        WriteI3Summary(summary, self.summaryfile)

        del tray
        return 0


   def Finish(self,stats={}):
        self.logger.info("finish %s: %s" % (self.__class__.__name__,
                                            self.GetParameter("execute")))
        return 0



class PPCResampleCorsika(PPC):
   """
    GPU Photon propagation w corsika resampling
   """

   def __init__(self):
        from I3Tray import I3Tray,I3Units
        PPC.__init__(self)
        self.AddParameter('OverSampling', 'Number of times to sample each shower.', 1)
        self.AddParameter('CylinderHeight', 'height of IceCube-centered target cylinder', 1200*I3Units.meter)
        self.AddParameter('CylinderRadius', 'radius of IceCube-centered target cylinder', 600*I3Units.meter)

   def Execute(self,stats={}):
        from .. import segments
        if not ipmodule.ParsingModule.Execute(self,stats): return 0

        from icecube import icetray, dataclasses, dataio, phys_services, interfaces
        from icecube import corsika_reader
        from I3Tray import I3Tray,I3Units

        # Instantiate a tray 
        tray = I3Tray()

        summary_in  = self.summaryfile
        summary_out = self.summaryfile
        if not os.path.exists(self.summaryfile):
           summary_in  = ''

        randomService = phys_services.I3SPRNGRandomService(
             		seed = self.seed,
             		nstreams = self.nproc,
             		streamnum = self.procnum)\
            if not self.usegslrng else phys_services.I3GSLRandomService(seed = self.seed*self.nproc+self.procnum)

        tray.context['I3RandomService'] = randomService



        inputfiles  = self.GetParameter('inputfilelist')
        # Configure IceTray modules 
        tray.AddModule("I3Reader", "reader",filenamelist=[self.gcdfile]+inputfiles)
        tray.AddModule("CORSIKAResampler","resample",
	        OverSampling=self.oversampling,
	        CylinderHeight=self.cylinderheight, 
	        CylinderRadius=self.cylinderradius)
        self.Configure(tray)

       # Execute 
        print(tray) 
        tray.Execute()

        del tray
        return 0




