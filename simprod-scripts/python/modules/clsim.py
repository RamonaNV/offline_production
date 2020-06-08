#!/usr/bin/env python
################################ Tray 1 ##############################
#   
#   IceProd classes for processing hit-generation in IceSim.
#   Author: juancarlos@icecube.wisc.edu
#   
######################################################################
import os
import sys,getopt,string
from os.path import expandvars
from .. import ipmodule
from ..util import CombineHits, DrivingTime, choose_max_efficiency
from ..util import ReadI3Summary, WriteI3Summary
CombineHits, DrivingTime, choose_max_efficiency
from ..util.fileutils import download,untar,isurl

def save_stats(tray,summary,stats):
    # Save stats
    for k in tray.Usage():
        stats[str(k.key())+":usr"] = k.data().usertime
        summary[str(k.key())+":usr"] = k.data().usertime

        stats[str(k.key())+":sys"] = k.data().systime
        summary[str(k.key())+":sys"] = k.data().systime

        stats[str(k.key())+":ncall"] = k.data().ncall
        summary[str(k.key())+":ncall"] = k.data().ncall



def LoadCascadeTables(Use_cvmfs=False, amplitude_fits = None, timing_fits = None):

     # Now we import the photon stuff
     from icecube import photonics_service

     # get the amplitude spline table
     if Use_cvmfs:
            spline_path = os.path.join(Use_cvmfs,'data','photon-tables','splines',
                                       os.path.basename(amplitude_fits))
     if Use_cvmfs and os.path.exists(spline_path):
            amplitude_fits = spline_path
     elif isurl(amplitude_fits):
            download(amplitude_fits)
            amplitude_fits = os.path.basename(amplitude_fits)
        
     # get the timing spline table
     if Use_cvmfs:
            spline_path = os.path.join(Use_cvmfs,'data','photon-tables','splines',
                                       os.path.basename(timing_fits))
     if Use_cvmfs and os.path.exists(spline_path):
            timing_fits = spline_path
     elif isurl(timing_fits):
            download(timing_fits)
            timing_fits = os.path.basename(timing_fits)

     return photonics_service.I3PhotoSplineService(
           amplitudetable=amplitude_fits,
           timingtable=timing_fits,
           timingSigma=0.)



####### IPModules ##################

class ClSim(ipmodule.ParsingModule):
   """
    GPU Photon propagation
   """
   def __init__(self):
        ipmodule.ParsingModule.__init__(self)

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
                          None)
        self.AddParameter("UseGPUs", "Use Graphics Processing Unit",False)
        self.AddParameter('RunMPHitFilter',"Run polyplopia's mphitfilter",True)
        self.AddParameter('PropagateMuons', 'Run PROPOSAL to do in-ice propagation', True)
        self.AddParameter('PROPOSALParams','any other parameters for proposal',dict())
        self.AddParameter("oversize","over-R: DOM radius oversize scaling factor",5)
        self.AddParameter("holeiceparametrization", "Location of hole ice param files", 
                          expandvars("$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm"))
        self.AddParameter("efficiency","overall DOM efficiency correction",[1.00])
        self.AddParameter("volumecyl","set volume to regular cylinder (set to False for 300m spacing from the DOMs)",True)
        self.AddParameter("IceModelLocation","Location of ice model param files", expandvars("$I3_BUILD/ice-models/resources/models")) 
        self.AddParameter("IceModel","ice model subdirectory", "spice_3.2") 
        self.AddParameter("MCTreeName","Name of MCTree frame object", "I3MCTree") 
        self.AddParameter("PhotonSeriesName","Photon Series Name","I3MCPESeriesMap") 
        self.AddParameter("RawPhotonSeriesName","Raw Photon Series Name",None) 
        self.AddParameter("UseGeant4","Enable Geant4 propagation",False) 
        self.AddParameter('HistogramFilename', 'Histogram filename.', None)
        self.AddParameter('EnableHistogram', 'Write a SanityChecker histogram file.', False)
        self.AddParameter('KeepMCTree','Delete propagated MCTree otherwise',True)
        self.AddParameter('KeepPureBackground','Keep events without signal hits',False)
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 

   def Configure(self,tray):

        from icecube import icetray, dataclasses, dataio, phys_services, interfaces
        from I3Tray import I3Tray,I3Units


        if self.gpu is not None and self.usegpus:
           os.putenv("CUDA_VISIBLE_DEVICES",str(self.gpu))
           os.putenv("COMPUTE",":0."+str(self.gpu))
           os.putenv("GPU_DEVICE_ORDINAL",str(self.gpu))

        # Configure IceTray services
        if os.path.exists(self.summaryfile):
           summary = ReadI3Summary(self.summaryfile)
        else:
           summary = dataclasses.I3MapStringDouble()
        tray.context['I3SummaryService'] = summary

        # RNG
        rngstate    = "rng.state"
        if not os.path.exists(rngstate):
           rngstate = ''
           self.logger.warning("Warning: no RNG state found. Using seed instead.")

        randomService = phys_services.I3SPRNGRandomService(
             		seed = self.seed,
             		nstreams = self.nproc,
             		streamnum = self.procnum)\
                        if not self.usegslrng else phys_services.I3GSLRandomService(seed = self.seed * self.nproc + self.procnum)
        tray.context['I3RandomService'] = randomService


        if type(self.efficiency) == list or type(self.efficiency) == tuple:
           if len(self.efficiency) == 1:
              efficiency=float(self.efficiency[0])
           elif len(self.efficiency) > 1:
              efficiency=map(float,self.efficiency)
           elif len(self.efficiency) > 1:
              raise Exception("Configured empty efficiency list")
        else:
            efficiency = choose_max_efficiency(self.efficiency)
        from .. import segments

        if self.propagatemuons:
        	randomServiceForPropagators = phys_services.I3SPRNGRandomService(
             		seed = self.seed,
             		nstreams = self.nproc*2,
             		streamnum = self.nproc + self.procnum)\
                        if not self.usegslrng else phys_services.I3GSLRandomService(seed = 2 * (self.seed * self.nproc + self.procnum))
        	tray.context['I3PropagatorRandomService'] = randomServiceForPropagators

        	tray.AddModule("Rename","rename_corsika_mctree",
        	           Keys=[self.mctreename,self.mctreename+'_preMuonProp'])
        	tray.AddSegment(segments.PropagateMuons, 'propagator',
                        RandomService= randomServiceForPropagators,
                        **self.proposalparams
        	) 

        from icecube import clsim
        tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
            GCDFile = self.gcdfile,
            RandomService = randomService,
            UseGPUs = self.usegpus,
            UseCPUs= not self.usegpus,
            IceModelLocation = os.path.join(self.icemodellocation,self.icemodel),
            UnshadowedFraction = efficiency,
            UseGeant4 = self.usegeant4,
            DOMOversizeFactor = self.oversize,
            MCTreeName = self.mctreename,
            MCPESeriesName = self.photonseriesname,
            PhotonSeriesName = self.rawphotonseriesname,
            HoleIceParameterization = self.holeiceparametrization
       	)


        if self.runmphitfilter:
            from icecube import polyplopia
            tray.AddModule("MPHitFilter","hitfilter",
                 HitOMThreshold=1,
                 RemoveBackgroundOnly= not self.keeppurebackground,
                 I3MCPESeriesMapName=self.photonseriesname)

        if self.enablehistogram and self.histogramfilename:         
            from icecube.production_histograms import ProductionHistogramModule
            from icecube.production_histograms.histogram_modules.simulation.mcpe_module import I3MCPEModule
        
            tray.AddModule(ProductionHistogramModule, 
                           Histograms = [I3MCPEModule],
                           OutputFilename = self.histogramfilename)

        if not self.keepmctree:
            self.logger.info("discarding %s" % (self.photonseriesname))
            tray.Add("Delete","clean_mctruth",
        	           Keys=[self.mctreename,self.mctreename+'_preSampling'])

        tray.AddModule("I3Writer","writer", 
                       filename=self.outputfile,
                       Streams=[icetray.I3Frame.TrayInfo,
                                icetray.I3Frame.DAQ,
                                icetray.I3Frame.Stream('S'),
                                icetray.I3Frame.Stream('M')])
         
   def Execute(self,stats):
        from .. import segments

        if not ipmodule.ParsingModule.Execute(self,stats): return 0

        from icecube import icetray, dataclasses, dataio, phys_services, interfaces
        from I3Tray import I3Tray,I3Units

        # Instantiate a tray 
        tray = I3Tray()

        inputfiles  = self.GetParameter('inputfilelist')
        # Configure IceTray modules 
        tray.AddModule("I3Reader", "reader",filenamelist=[self.gcdfile]+inputfiles)
        self.Configure(tray)

       # Execute the Tray
        tray.Execute()

        summary = tray.context['I3SummaryService']
        save_stats(tray,summary,stats)
        WriteI3Summary(summary, self.summaryfile)

        del tray
        return 0


   def Finish(self,stats={}):
        self.logger.info("finish %s: %s" % (self.__class__.__name__,self.GetParameter("execute")))
        return 0



class ClSimResampleCorsika(ClSim):
   """
    GPU Photon propagation w corsika resampling
   """

   def __init__(self):
        from I3Tray import I3Tray,I3Units
        ClSim.__init__(self)
        self.AddParameter('OverSampling', 'Number of times to sample each shower.', 1)
        self.AddParameter('CylinderHeight', 'height of IceCube-centered target cylinder', 1200*I3Units.meter)
        self.AddParameter('CylinderRadius', 'radius of IceCube-centered target cylinder', 600*I3Units.meter)

   def Execute(self,stats):
        from .. import segments
        if not ipmodule.ParsingModule.Execute(self,stats): return 0

        from icecube import icetray, dataclasses, dataio, phys_services, interfaces
        from icecube import corsika_reader
        from I3Tray import I3Tray,I3Units

        # Instantiate a tray 
        tray = I3Tray()

        inputfiles  = self.GetParameter('inputfilelist')
        # Configure IceTray modules 
        tray.AddModule("I3Reader", "reader",filenamelist=[self.gcdfile]+inputfiles)
        tray.AddModule("CORSIKAResampler","resample",
	        OverSampling=self.oversampling,
	        CylinderHeight=self.cylinderheight, 
	        CylinderRadius=self.cylinderradius)
        self.Configure(tray)

       # Execute the Tray
        tray.Execute()
        summary = tray.context['I3SummaryService']
        save_stats(tray,summary,stats)
        WriteI3Summary(summary, self.summaryfile)

        del tray
        return 0




class HybridPhotons(ipmodule.ParsingModule):
   """
    GPU+Spline-based Photon propagation
   """
   def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('Hybrid','Run in hybrid mode',True)
        self.AddParameter('KeepMCTree','Delete propagated MCTree otherwise',True)
        self.AddParameter('PropagateMuons', 'Run PROPOSAL to do in-ice propagation', True)
        self.AddParameter('PROPOSALParams','any other parameters for proposal',dict())
        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('inputfilelist','list of input filenames',[])
        self.AddParameter('outputfile','Output filename','')
        self.AddParameter('seed','RNG seed',0)
        self.AddParameter('procnum','job number (RNG stream number)',0)
        self.AddParameter('nproc','Number of jobs (Number of RNG streams)',1)
        self.AddParameter('summaryfile','JSON Summary filename','$steering(summaryfile)')
        self.AddParameter('MJD','MJD (0=do not modify)',0)
        self.AddParameter("GPU", 
                          "Graphics Processing Unit number (shoud default to environment if None)",
                          None)
        self.AddParameter("UseGPUs", "Use Graphics Processing Unit",False)
        self.AddParameter('RunMPHitFilter',"Run polyplopia's mphitfilter",True)
        self.AddParameter("oversize","over-R: DOM radius oversize scaling factor",5)
        self.AddParameter("efficiency","overall DOM efficiency correction",1.00)
        self.AddParameter("volumecyl","set volume to regular cylinder (set to False for 300m spacing from the DOMs)",True)
        self.AddParameter("IceModelLocation","Location of ice model param files", expandvars("$I3_BUILD/clsim/resources/ice")) 
        self.AddParameter("IceModel","ice model subdirectory", "SpiceMie") 
        self.AddParameter("PhotonSeriesName","Photon Series Name","I3MCPESeriesMap") 
        self.AddParameter("CscdSplineAmplitudeFits",
                        "Cacade spline amplitude .fits tables",
                        'ems_spice1_z20_a10.abs.fits')
        self.AddParameter("CscdSplineTimingFits",
                        "Cacade spline timing .fits tables",
                        'ems_spice1_z20_a10.prob.fits')
        self.AddParameter("usecvmfs","Use CVMFS for spline tables (if possible)",
                          '/cvmfs/icecube.opensciencegrid.org')
        self.AddParameter("ConvertToMCHits", "Convert MCPEs to MCHits for backwards compatibility", False)
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 

   def Execute(self,stats):
        if not ipmodule.ParsingModule.Execute(self,stats): return 0
        from ..util.fileutils import download,untar,isurl
        from icecube import icetray, dataclasses, dataio, phys_services, interfaces, sim_services
        from I3Tray import I3Tray, I3Units

        # Instantiate a tray 
        tray = I3Tray()

        gcdfile     = self.GetParameter('gcdfile')
        inputfiles  = self.GetParameter('inputfilelist')
        outputfile  = self.GetParameter('outputfile')


        # Configure IceTray services
        if os.path.exists(self.summaryfile):
           summary = ReadI3Summary(self.summaryfile)
        else: 
           summary = dataclasses.I3MapStringDouble()
        tray.context['I3SummaryService'] = summary

        # RNG
        rngstate    = "rng.state"
        if not os.path.exists(rngstate): 
           rngstate = ''
           self.logger.warning("no RNG state found. Using seed instead.")

        randomService = phys_services.I3SPRNGRandomService(
             		seed = self.seed,
             		nstreams = self.nproc,
             		streamnum = self.procnum)\
            if not self.usegslrng else phys_services.I3GSLRandomService(seed = self.seed * self.nproc + self.procnum)

        tray.context['I3RandomService'] = randomService


        # Configure IceTray modules 
        tray.AddModule("I3Reader", "reader",filenamelist=[gcdfile]+inputfiles)
        
        # simulate cascades (with photonics-hit-maker)
        cascade_service = LoadCascadeTables(
            Use_cvmfs=self.usecvmfs, 
            amplitude_fits = self.cscdsplineamplitudefits, 
            timing_fits = self.cscdsplinetimingfits)

        from .. import segments
        if self.propagatemuons:
        	randomServiceForPropagators = phys_services.I3SPRNGRandomService(
             		seed = self.seed,
             		nstreams = self.nproc*2,
             		streamnum = self.nproc + self.procnum)\
                        if not self.usegslrng else phys_services.I3GSLRandomService(seed = 2 * (self.seed * self.nproc + self.procnum))
        	tray.context['I3PropagatorRandomService'] = randomServiceForPropagators

        	tray.AddModule("Rename","rename_corsika_mctree",Keys=['I3MCTree','I3MCTree_preMuonProp'])
        	tray.AddSegment(segments.PropagateMuons, 'propagator',
                        RandomService= randomServiceForPropagators,
                        **self.proposalparams
        	) 

        tray.AddSegment(segments.PropagatePhotons, "hybridpes",
            GCDFile=gcdfile,
            RandomService = randomService,
            HybridMode = True,
            UseGPUs = False,
            UseAllCPUCores = False,
            IceModel = self.icemodel,
            IceModelLocation = self.icemodellocation,
            CascadeService = cascade_service,
            UseCascadeExtension = False,
            UnshadowedFraction = self.efficiency,
            DOMOversizeFactor = self.oversize,
            OutputPESeriesMapName = self.photonseriesname)


        if self.runmphitfilter:
           from icecube import polyplopia
           tray.AddModule("MPHitFilter","hitfilter",
                 HitOMThreshold=1,
                 RemoveBackgroundOnly=False,
                 I3MCPESeriesMapName=self.photonseriesname)

        # Convert for backwards compatibility
        if self.GetParameter("ConvertToMCHits"):
           tray.AddModule("I3MCPEtoI3MCHitConverter","converthits", 
                          InputResponse=self.photonseriesname,
                          OutputResponse = "MCHitSeriesMap")
           tray.AddModule( DrivingTime, "dt", Streams = [icetray.I3Frame.DAQ] )

        if not self.keepmctree:
            self.logger.info("discarding %s" % (self.photonseriesname))
            tray.Add("Delete","clean_mctruth",
        	           Keys=[self.mctreename,self.mctreename+'_preSampling'])

        tray.AddModule("I3Writer","writer", filename=outputfile,streams=[icetray.I3Frame.DAQ])
         
        # Execute the Tray
        tray.Execute()

        summary = tray.context['I3SummaryService']
        save_stats(tray,summary,stats)
        WriteI3Summary(summary, self.summaryfile)

        del tray
        return 0


   def Finish(self,stats={}):
        self.logger.info("finish %s: %s" % (self.__class__.__name__,self.GetParameter("execute")))
        return 0


