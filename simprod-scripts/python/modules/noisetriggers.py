#!/usr/bin/env python
##################################################################################
# A script designed to produce events containing only noise triggers. Contains
# both a normal tray segment (with everything except a random module and an
# I3Writer) and an iceprod module for SimProd generation.
#
# Each event produced will initally be a 100 ms long frame of noise. This will
# be triggered and then cut down using CoincidenceAfterProcessing from trigger-sim
# to create normal-sized events.
#
# @author Michael Larson (mjlarson@nbi.ku.dk)
# @date 15 July 2015
##################################################################################
from ..segments import ProduceNoiseTriggers
from ..util import ReadI3Summary, WriteI3Summary
from .. import ipmodule
from I3Tray import I3Tray
from icecube import icetray, dataio, phys_services, dataclasses
import os


class NoiseTriggers(ipmodule.ParsingModule):
    def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('summaryfile','JSON Summary filename','summary.json')
        self.AddParameter('outputfile','Output filename','')
        self.AddParameter('mjd','MJD for the GCD file',56063)
        self.AddParameter("Detector","detector to simulate","IC86:2011")
        self.AddParameter('nevents','Number of events',1000)
        self.AddParameter("RunId","Configure run ID",0)
        self.AddParameter("RNGSeed","RNG seed",0)

    def Execute(self, stats):
        if not ipmodule.ParsingModule.Execute(self, stats): return 0

        tray = I3Tray()
        randomService = phys_services.I3GSLRandomService(self.rngseed)
        tray.context['I3RandomService'] = randomService

        summary = dataclasses.I3MapStringDouble()
        summaryfile     = self.GetParameter('summaryfile')
        if os.path.exists(summaryfile): 
           summary = ReadI3Summary(summaryfile)
        tray.context['I3SummaryService'] = summary

        tray.AddSegment(ProduceNoiseTriggers, "noise_triggers",
                        gcd_file = self.gcdfile,
                        nevents = self.nevents,
                        run_id = self.runid)

        from ..util import BasicCounter, DAQCounter
        from icecube import icetray
        tray.AddModule(BasicCounter,"count_triggers",
                       Streams = [icetray.I3Frame.DAQ] ,
                       name="%s Triggered Events" % self.GetParameter('Detector'),
                       Stats=stats)
        
        skipkeys = [ "I3Triggers", "MCPMTResponseMap", "MCTimeIncEventID"]

        tray.AddModule("I3Writer","writer",
                       filename=self.outputfile,
                       Streams=[icetray.I3Frame.TrayInfo,
                                icetray.I3Frame.DAQ,
                                icetray.I3Frame.Stream('S'),
                                icetray.I3Frame.Stream('M')],
                       SkipKeys = skipkeys,
                       )

        tray.Execute(4+self.nevents)
        
        tray.PrintUsage()

        summary = tray.context['I3SummaryService']
        # Save stats
        for k in tray.Usage():
            stats[str(k.key())+":usr"] = k.data().usertime
            summary[str(k.key())+":usr"] = k.data().usertime

            stats[str(k.key())+":sys"] = k.data().systime
            summary[str(k.key())+":sys"] = k.data().systime

            stats[str(k.key())+":ncall"] = k.data().ncall
            summary[str(k.key())+":ncall"] = k.data().ncall

        WriteI3Summary(summary, summaryfile)

        return 0

