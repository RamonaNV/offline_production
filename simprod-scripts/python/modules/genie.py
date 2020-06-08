#!/usr/bin/env python

# This a working script taken from Ken for generating PINGU.
# The goal is to create something that can run using the latest
# release of simulation instead of the non-java decaf version


from I3Tray import *
import os
import os.path
import sys
import random

from .. import ipmodule
from ..util import ReadI3Summary, WriteI3Summary


def AddEmptyWeights(frame): 
   frame["WeightDict"] = dataclasses.I3MapStringDouble()
   
def ResetMCPE(frame):
  frame.Delete("I3MCPESeriesMap")
  frame["I3MCPESeriesMap"] = frame["CleanedI3MCPESeriesMap"]
  frame.Delete("CleanedI3MCPESeriesMap")

# Genie-icetray doesn't respect the absolute position of the GenVolDepth. Correct the x/y/z position now.
def shifter(frame, x=0, y=0, z=0):
   from icecube import dataclasses
   new_center = dataclasses.I3Position(x * I3Units.m, y * I3Units.m, z * I3Units.m)
   mctree = frame["I3MCTree"]
   particle = mctree.get_head()
   
   # Recursively shift all of the particles in the MC tree to the new coordinates
   def shift_daughters(tree, p):
      daughters = tree.get_daughters(p)
      if len(daughters) > 0:
         for d in daughters: shift_daughters(tree, d)
         tree[p].pos += new_center
         return

   shift_daughters(mctree, mctree.get_head())
   del frame["I3MCTree"]
   frame["I3MCTree"] = mctree
   return


class Genie(ipmodule.ParsingModule):

   def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('geniepath',"GENIE's very own ROOTSYS",'')
        self.AddParameter('outputfile','Output filename','')
        self.AddParameter('summaryfile','JSON Summary filename','genie.json')
        self.AddParameter('nevents','Number of events',1000)
        self.AddParameter('emax','Maximum energy',190.0*I3Units.GeV)
        self.AddParameter('emin','Minimun energy',1.0*I3Units.GeV)
        self.AddParameter('gamma','Gamma index',1.0)
        self.AddParameter('NuFlavor','Neutrino Flavor','NuE')
        self.AddParameter("seed","RNG seed",0)
        self.AddParameter("procnum","RNG stream number",0)
        self.AddParameter("nproc","Number of RNG streams",1)
        self.AddParameter('Polyplopia','Produce coincident showers',False)
        self.AddParameter('BackgroundFile','pre-generated coincident showers file',"")
        self.AddParameter("length", "cylinder length in m", 1200.)
        self.AddParameter("radius", "cylinder radius in m", 800.)
        self.AddParameter("x", "cylinder x-position in m", 0.)
        self.AddParameter("y", "cylinder y-position in m", 0.)
        self.AddParameter("z", "cylinder z-position in m", 0.)
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 

   def Execute(self,stats):
        if not ipmodule.ParsingModule.Execute(self,stats): return 0

        from icecube import icetray
        from icecube import dataclasses
        from icecube import dataio
        from icecube import phys_services
        from icecube import sim_services
        from icecube import genie_icetray
        from icecube import simclasses


        # Instantiate a tray 
        tray = I3Tray()

        # Configure IceTray services
        summary = dataclasses.I3MapStringDouble()
        tray.context['I3SummaryService'] = summary

        randomService = phys_services.I3SPRNGRandomService(
             seed = self.seed, 
             nstreams = self.nproc, 
             streamnum = self.procnum)\
            if not self.usegslrng else phys_services.I3GSLRandomService(seed = self.seed*self.nproc+self.procnum)

        randomServiceForPropagators = phys_services.I3SPRNGRandomService(
             seed = self.seed,
             nstreams = self.nproc*2,
             streamnum = self.nproc+ self.procnum)\
            if not self.usegslrng else phys_services.I3GSLRandomService(seed = 2*(self.seed*self.nproc+self.procnum))

        tray.context['I3RandomService'] = randomService
        tray.context['I3PropagatorRandomService'] = randomServiceForPropagators

        tray.AddModule("I3InfiniteSource","TheSource",
                       Prefix=self.gcdfile,
                       Stream=icetray.I3Frame.DAQ)
    
        tray.AddModule("I3GENIEGenerator","genie_generator",
 	       RandomService   = randomService, # alternatively, this can be None and the I3RandomService can be installed using tray.AddService()
 	       SplineFilename  = os.path.expandvars("$I3_BUILD/genie-icetray/resources/splines/splines_water_2.8.6.xml"),
 	       LHAPDFPath      = os.path.expandvars("$I3_BUILD/genie-icetray/resources/PDFsets"), # use $I3_PORTS/bin/lhapdf-getdata to download the PDFsets
 	       NuEnergyMin          = self.emin,
 	       NuEnergyMax          = self.emax,
 	       PowerLawIndex        = self.gamma, 
 	       GenVolRadius         = self.radius*I3Units.m,
 	       GenVolLength         = self.length*I3Units.m,
 	       GenVolDepth          = 1950.*I3Units.m,
	       NeutrinoFlavor       = self.nuflavor, # generates neutrinos and anti-neutrinos (1:1)
 	       MaterialDensity      = 0.93*I3Units.g/I3Units.cm3, # ice density
 	       TargetMixIngredients = [1000080160,1000010010], # O16, H1
	       TargetMixQuantities  = [1,2], # H2O (O16->1x, H1->2x)
 	       ForceSingleProbScale = False,
 	       NEvents              = self.nevents,
 	       GENIEPATH            = self.geniepath
 	)

        tray.AddModule(shifter, "cylinder_shifter",
                       x = self.x * I3Units.m,
                       y = self.y * I3Units.m,
                       z = self.z * I3Units.m,
                       Streams=[icetray.I3Frame.DAQ])           

        # Add empty weightdict 
        tray.AddModule(AddEmptyWeights)

        if self.polyplopia:

            from .. import segments
            tray.AddSegment(segments.PolyplopiaSegment,"coincify",
                    RandomService = randomService,
                    mctype = 'Genie',
                    mctree_name = "I3MCTree",
                    separate_coincident_mctree_name = "CoincidentI3MCTree_preMuonProp",
                    bgfile = self.backgroundfile,
                    timewindow = 40.*I3Units.microsecond,
                    rate = 5.0*I3Units.kilohertz,
            ) 

            tray.AddSegment(segments.PropagateMuons, 'propagator',
                    RandomService = randomServiceForPropagators,
                    InputMCTreeName ="CoincidentI3MCTree_preMuonProp",
                    OutputMCTreeName ="CoincidentI3MCTree",
            ) 
 
        tray.AddModule("I3Writer","writer",
            filename = self.outputfile,
            Streams=[icetray.I3Frame.TrayInfo,
                     icetray.I3Frame.DAQ,
                     icetray.I3Frame.Stream('S'),
                     icetray.I3Frame.Stream('M')])

        
        tray.Execute()
        
        return 0



class GeniePlusClSim(ipmodule.ParsingModule):

   def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('geniepath',"GENIE's very own ROOTSYS",'')
        self.AddParameter('outputfile','Output filename','')
        self.AddParameter('summaryfile','Summary filename','genie.json')
        self.AddParameter('nevents','Number of events',1000)
        self.AddParameter('gpuopt','GPU',0)
        self.AddParameter('emax','Maximum energy',190.0*I3Units.GeV)
        self.AddParameter('emin','Minimun energy',1.0*I3Units.GeV)
        self.AddParameter('gamma','Gamma index',1.0)
        self.AddParameter('NuFlavor','Neutrino Flavor','NuE')
        self.AddParameter("RNGSeed","RNG seed",0)
        self.AddParameter("RNGStream","RNG stream number",0)
        self.AddParameter("RNGNumberOfStreams","Number of RNG streams",1)

        self.AddParameter("length", "cylinder length in m", 1200.)
        self.AddParameter("radius", "cylinder radius in m", 800.)
        self.AddParameter("x", "cylinder x-position in m", 0.)
        self.AddParameter("y", "cylinder y-position in m", 0.)
        self.AddParameter("z", "cylinder z-position in m", 0.)

        self.AddParameter("oversize","over-R: DOM radius oversize scaling factor",1)
        self.AddParameter("efficiency","overall DOM efficiency correction",0.99)
        self.AddParameter("IceModelLocation","Location of ice model param files", os.path.expandvars("$I3_BUILD/clsim/resources/ice")) 
        self.AddParameter("IceModel","ice model subdirectory", "SpiceMie") 
        self.AddParameter("PhotonSeriesName","Photon Series Name","I3MCPESeriesMap") 
        self.AddParameter("RawPhotonSeriesName","Raw Photon Series Name",None)
        self.AddParameter('Polyplopia','Produce coincident showers',False)
        self.AddParameter('BackgroundFile','pre-generated coincident showers file',"")
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 

   def Execute(self,stats):
        if not ipmodule.ParsingModule.Execute(self,stats): return 0

        from icecube import icetray
        from icecube import dataclasses
        from icecube import dataio
        from icecube import phys_services
        from icecube import sim_services
        from icecube import genie_icetray
        from icecube import simclasses
        from ..segments import PropagatePhotons
        from ..util import BasicHitFilter
        from I3Tray import I3Tray,I3Units


        icetray.set_log_level_for_unit('I3CLSimStepToPhotonConverterOpenCL', icetray.I3LogLevel.LOG_FATAL)

        # Instantiate a tray 
        tray = I3Tray()

        randomService = phys_services.I3SPRNGRandomService(self.rngseed, self.rngnumberofstreams, self.rngstream)\
            if not self.usegslrng else phys_services.I3GSLRandomService(seed = self.rngseed*self.rngnumberofstreams+self.rngstream)

        randomServiceForPropagators = phys_services.I3SPRNGRandomService(
             seed = self.rngseed,
             nstreams = self.rngnumberofstreams*2,
             streamnum = self.rngnumberofstreams + self.rngstream)\
            if not self.usegslrng else phys_services.I3GSLRandomService(seed = 2*(self.rngseed*self.rngnumberofstreams+self.rngstream))

        tray.context['I3RandomService'] = randomService
        tray.context['I3PropagatorRandomService'] = randomServiceForPropagators

        # Configure IceTray services
        summary = dataclasses.I3MapStringDouble()
        tray.context['I3SummaryService'] = summary


        tray.AddModule("I3InfiniteSource","TheSource",
                       Prefix=self.gcdfile,
                       Stream=icetray.I3Frame.DAQ)
    
        tray.AddModule("I3GENIEGenerator","genie_generator",
 	       RandomService   = randomService, # alternatively, this can be None and the I3RandomService can be installed using tray.AddService()
 	       SplineFilename  = os.path.expandvars("$I3_BUILD/genie-icetray/resources/splines/splines_water_2.8.6.xml"),
 	       LHAPDFPath      = os.path.expandvars("$I3_BUILD/genie-icetray/resources/PDFsets"), # use $I3_PORTS/bin/lhapdf-getdata to download the PDFsets
 	       NuEnergyMin          = self.emin,
 	       NuEnergyMax          = self.emax,
 	       PowerLawIndex        = self.gamma, 
 	       GenVolRadius         = self.radius*I3Units.m,
 	       GenVolLength         = self.length*I3Units.m,
 	       GenVolDepth          = 1950.*I3Units.m,
	       NeutrinoFlavor       = self.nuflavor, # generates neutrinos and anti-neutrinos (1:1)
 	       MaterialDensity      = 0.93*I3Units.g/I3Units.cm3, # ice density
 	       TargetMixIngredients = [1000080160,1000010010], # O16, H1
	       TargetMixQuantities  = [1,2], # H2O (O16->1x, H1->2x)
 	       ForceSingleProbScale = False,
 	       NEvents              = self.nevents,
 	       GENIEPATH            = self.geniepath
 	)

        tray.AddModule(shifter, "cylinder_shifter",
                       x = self.x * I3Units.m,
                       y = self.y * I3Units.m,
                       z = self.z * I3Units.m,
                       Streams=[icetray.I3Frame.DAQ])           

        tray.AddSegment(PropagatePhotons,"photons",
			UseGPUs = self.gpuopt,
			UnshadowedFraction = self.efficiency,
			DOMOversizeFactor = self.oversize,
			IceModelLocation = self.icemodellocation,
                        HybridMode = False,
			IceModel = self.icemodel,
			OutputPESeriesMapName = self.photonseriesname,
			OutputPhotonSeriesName = self.rawphotonseriesname,
			UseGeant4 = True,
		)
        tray.AddModule(BasicHitFilter,"hitfilter",
			Streams = [icetray.I3Frame.DAQ],
			HitSeriesMapNames = self.photonseriesname)

        # Add empty weightdict 
        tray.AddModule(AddEmptyWeights)
        tray.AddModule("I3RemoveLargeDT","removeLDT")
        tray.AddModule(ResetMCPE,'reset',Streams=[icetray.I3Frame.DAQ])

        if self.polyplopia:

            from .. import segments
            tray.AddSegment(segments.PolyplopiaSegment,"coincify",
                    RandomService=randomService,
                    mctype='Genie',
                    mctree_name="I3MCTree",
                    separate_coincident_mctree_name="CoincidentI3MCTree_preMuonProp",
                    bgfile = self.backgroundfile,
                    timewindow = 40.*I3Units.microsecond,
                    rate = 5.0*I3Units.kilohertz,
            ) 

            tray.AddSegment(segments.PropagateMuons, 'propagator',
                    RandomService = randomServiceForPropagators,
                    InputMCTreeName = "CoincidentI3MCTree_preMuonProp",
                    OutputMCTreeName = "CoincidentI3MCTree",
            ) 

        tray.AddModule("I3Writer","writer",
            filename = self.outputfile,
            Streams=[icetray.I3Frame.TrayInfo,
                     icetray.I3Frame.DAQ,
                     icetray.I3Frame.Stream('S'),
                     icetray.I3Frame.Stream('M')])


        
        tray.Execute()
        
        return 0
