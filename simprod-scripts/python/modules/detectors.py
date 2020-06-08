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
from ..util import ReadI3Summary, WriteI3Summary

logger = logging.getLogger('detector')


class IceCube(ipmodule.ParsingModule):

    def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('SummaryFile','Summary filename', None)
        self.AddParameter('MCType','Generator particle type','corsika_weighted')
        self.AddParameter('UseLinearTree','Use I3LinearizedMCTree for serialization', False)
        self.AddParameter('MCPrescale','Prescale for keeping additional Monte Carlo info in the frame', 100)
        self.AddParameter('IceTop','Do IceTop Simulation?', False)
        self.AddParameter('Genie','Assume separate Genie MCPEs and BG MCPEs', False)
        self.AddParameter('FilterTrigger','filter untriggered events', True)
        self.AddParameter('Trigger','Run trigger simulation', True)
        self.AddParameter('LowMem','Low Memory mode', False)
        self.AddParameter('BeaconLaunches','Simulate beacon launches', True)
        self.AddParameter('TimeShiftSkipKeys','Skip keys in the triggersim TimeShifter', [])
        self.AddParameter('SampleEfficiency','Resample I3MCPESeriesMap for different efficiency', 0.0)
        self.AddParameter('GeneratedEfficiency','Generated efficiency for resampling', 0.0)
        self.AddParameter('RunID','Run ID', 0, explicit_type='int')
        self.AddParameter('gcdfile', 'GeoCalibDetStatus filename', '')
        self.AddParameter('inputfile', 'Input filename', '')
        self.AddParameter('seed', 'RNG Seed', 123)
        self.AddParameter('procnum', 'RNG stream number', 0)
        self.AddParameter('nproc', 'RNG number of streams', 1)
        self.AddParameter('MCPESeriesName', 'Name of MCPESeriesMap in frame', 'I3MCPESeriesMap')
        self.AddParameter('HistogramFilename', 'Histogram filename.', None)
        self.AddParameter('EnableHistogram', 'Write a SanityChecker histogram file.', False)
        self.AddParameter('outputfile', 'Output filename', 'output.i3.gz')
        self.AddParameter('DetectorName', 'Name of detector', 'IC86')
        self.AddParameter('SkipKeys', 'Skip keys for the writer', [])
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 




    def Execute(self, stats):
        'Runs the tray, but the important stuff happens in segment_X(...)'
        if not ipmodule.ParsingModule.Execute(self,stats):
            return 0

        # Load libraries
        from icecube import clsim
        from icecube import phys_services
        from icecube import sim_services
        from icecube import vuvuzela
        from icecube import DOMLauncher
        from icecube import trigger_sim

        # Instantiate a tray
        tray = I3Tray()

        randomService = phys_services.I3SPRNGRandomService(
            seed = self.seed,
            nstreams = self.nproc,
            streamnum = self.procnum)\
            if not self.usegslrng else phys_services.I3GSLRandomService(seed = self.seed*self.nproc+self.procnum)

        tray.context['I3RandomService'] = randomService

        # Configure IceTray modules
        tray.AddModule("I3Reader","reader", FilenameList=[self.gcdfile, self.inputfile])

        if self.mcpeseriesname != 'I3MCPESeriesMap':
            def mover(fr):
                if 'I3MCPESeriesMap' in fr:
                    del fr['I3MCPESeriesMap']
                fr['I3MCPESeriesMap'] = fr[self.mcpeseriesname]
            tray.Add(mover, "move_MCPESeries",Streams=[icetray.I3Frame.DAQ])


        # Instantiate a SummaryService if required
        summary = dataclasses.I3MapStringDouble()
        if self.summaryfile and os.path.exists(self.summaryfile):
              summary = ReadI3Summary(self.summaryfile)
        tray.context['I3SummaryService'] = summary

        
        tray.AddSegment(segments.DetectorSegment,"detector",
            gcdfile=self.gcdfile,
            mctype=self.mctype,
            uselineartree=self.uselineartree,
            detector_label=self.detectorname,
            runtrigger=self.trigger,
            filtertrigger=self.filtertrigger,
            stats=stats,
            icetop=self.icetop,
            genie=self.genie,
            prescale=self.mcprescale,
            lowmem=self.lowmem,
            BeaconLaunches=self.beaconlaunches,
            TimeShiftSkipKeys=self.timeshiftskipkeys,
            SampleEfficiency=self.sampleefficiency,
            GeneratedEfficiency=self.generatedefficiency,
            RunID=self.runid,
            KeepMCHits = not self.procnum % self.mcprescale,
            KeepPropagatedMCTree = not self.procnum % self.mcprescale,
            KeepMCPulses = not self.procnum % self.mcprescale
        )

        if self.enablehistogram and self.histogramfilename:         
            from icecube.production_histograms import ProductionHistogramModule
            from icecube.production_histograms.histogram_modules.simulation.pmt_response import PMTResponseModule
            from icecube.production_histograms.histogram_modules.simulation.dom_mainboard_response import InIceResponseModule
            from icecube.production_histograms.histogram_modules.simulation.trigger import TriggerModule
            from icecube.production_histograms.histograms.simulation.noise_occupancy import NoiseOccupancy

            tray.AddModule(ProductionHistogramModule, 
                           Histograms = [PMTResponseModule,
                                         InIceResponseModule,
                                         TriggerModule,
                                         NoiseOccupancy],
                           OutputFilename = self.histogramfilename)

        tray.AddModule("I3Writer","writer",
            Filename=self.outputfile,
            SkipKeys=self.skipkeys,
            Streams=[icetray.I3Frame.TrayInfo,
                     icetray.I3Frame.DAQ,
                     icetray.I3Frame.Stream('S'),
                     icetray.I3Frame.Stream('M')])



        # Execute the Tray
        tray.Execute()

        for k in tray.Usage():
            stats[str(k.key())+":usr"] = k.data().usertime
            stats[str(k.key())+":sys"] = k.data().systime
            stats[str(k.key())+":ncall"] = k.data().systime

        # Free memory
        del tray
        return 0



class IceTop(ipmodule.ParsingModule):
  """
  Iceprod module to wrap iceprod.segments.IceTopSim. This is a convenience
  module that makes a simple simulation of the detector for IceTop only.
  """
  def __init__(self):
    ipmodule.ParsingModule.__init__(self)

    self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
    self.AddParameter('inputfile','Input filename','')
    self.AddParameter('outputfile','Output filename','')
    self.AddParameter('SummaryFile','Summary filename', None)
    self.AddParameter('Seed','RNG Seed', 123)
    self.AddParameter('procnum','RNG stream number', 0)
    self.AddParameter('nproc','RNG number of streams', 1)
    self.AddParameter('DOMLauncher', 'Simulate with DOMLauncher', True)
    self.AddParameter('sim_trigger', 'Simulate trigger', False)
    self.AddParameter('calibrate', 'Calibrate and extract pulses (requires tpx module, which is in IceRec usually)', False)
    self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 

  def Execute(self,stats):
    if not ipmodule.ParsingModule.Execute(self,stats):
      return 0

    from icecube import phys_services
    from icecube import sim_services
    from .. import segments

    tray = I3Tray()

    summary = dataclasses.I3MapStringDouble()
    if self.summaryfile and os.path.exists(self.summaryfile):
       summary = ReadI3Summary(self.summaryfile)
    tray.context['I3SummaryService'] = summary


    # set the random number generator
    rngstate    = "rng.state"
    if not os.path.exists(rngstate):
      rngstate = ''
      self.logger.warning("no RNG state found. Using seed instead.")

    if not self.usegslrng:
        tray.AddService("I3SPRNGRandomServiceFactory","random",
                        Seed = self.seed,
                        StreamNum = self.procnum,
                        NStreams = self.nproc,
                        instatefile = rngstate,
                        outstatefile = 'rng.state'
                    )
    else:    
        tray.context["I3RandomService"] = phys_services.I3GSLRandomService(seed)


    tray.AddModule("I3Reader","reader",
                   filenamelist = [self.gcdfile, self.inputfile],
                   )

    # The main segment
    tray.AddSegment(segments.IceTopSim, 'IceTopSim',
                    sim_trigger = self.sim_trigger,
                    sim_new = self.domlauncher,
                    calibrate = self.calibrate,
                    gcd = self.gcdfile,
                    )

    tray.AddModule("I3Writer","writer",
                   filename = self.outputfile,
                   Streams=[icetray.I3Frame.TrayInfo,
                            icetray.I3Frame.DAQ,
                            icetray.I3Frame.Stream('S'),
                            icetray.I3Frame.Stream('M')],
                   )

    tray.Execute()

    summary = tray.context['I3SummaryService']
    for k in tray.Usage():
            stats[str(k.key())+":usr"] = k.data().usertime
            summary[str(k.key())+":usr"] = k.data().usertime

            stats[str(k.key())+":sys"] = k.data().systime
            summary[str(k.key())+":sys"] = k.data().systime

            stats[str(k.key())+":ncall"] = k.data().ncall
            summary[str(k.key())+":ncall"] = k.data().ncall

    if self.summaryfile:
       WriteI3Summary(summary, self.summaryfile)

    # Free memory
    del tray
    return 0

