from os.path import expandvars
from icecube.simprod import ipmodule


class MuonPropagator(ipmodule.ParsingModule):
    """IceProd module for ``MuonGun`` simulations

    It defines parameters for and executes
    :func:`icecube.simprod.segments.GenerateCosmicRayMuons` and
    :func:`icecube.simprod.segments.PropagateMuons`.

    **Module parameters:**

        :attr:`inputfilelist`
            List of I3 files to read and propagate
        :attr:`nproc`
            Number of processes for RNG
        :attr:`pronum`
            Process number
        :attr:`seed`
            RNG seed
        :attr:`length`
            Cylinder length in m (now depricated, use config_file)
        :attr:`radius`
            Cylinder radius in m (now depricated, use config_file)
        :attr:`x`
            Cylinder x-position in m (now depricated, use config_file)
        :attr:`y`
            Cylinder y-position in m (now depricated, use config_file)
        :attr:`z`
            Cylinder z-position in m (now depricated, use config_file)
        :attr:`proposalconfigfile`
            path to configuration file for PROPOSAL
        :attr:`outputfile`
            Output filename
        :attr:`summaryfile`
            Summary filename

    """
    def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter("nproc", "number of processes for RNG", 1)
        self.AddParameter("procnum", "process number", 0)
        self.AddParameter("seed", "RNG seed", 1)

        self.AddParameter("length", "cylinder length in m (now depricated, use config_file)", 1600.)
        self.AddParameter("radius", "cylinder radius in m (now depricated, use config_file)", 800.)
        self.AddParameter("x", "cylinder x-position in m (now depricated, use config_file)", 0.)
        self.AddParameter("y", "cylinder y-position in m (now depricated, use config_file)", 0.)
        self.AddParameter("z", "cylinder z-position in m (now depricated, use config_file)", 0.)
        self.AddParameter("z", "cylinder z-position in m (now depricated, use config_file)", 0.)
        self.AddParameter("proposalconfigfile",
                          "path to configuration file for PROPOSAL",
                          expandvars("$I3_BUILD/PROPOSAL/resources/config_icesim.json"))

        self.AddParameter('inputfilelist','list of input filenames',[])
        self.AddParameter('InputMCTreeName','name of MCTree object in frame',"I3MCTree_no_propagation",)
        self.AddParameter('OutputMCTreeName','name of MCTree object to be added to frame',"I3MCTree")
        self.AddParameter("outputfile", "output filename", "corsika.i3")
        self.AddParameter('summaryfile','JSON Summary filename','summary.json')

        self.AddParameter('HistogramFilename', 'Histogram filename.', None)
        self.AddParameter('EnableHistogram', 'Write a SanityChecker histogram file.', False)
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 

    def Execute(self, stats):
        if not ipmodule.ParsingModule.Execute(self, stats):
            return 0

        import icecube.icetray
        import icecube.dataclasses
        import icecube.dataio
        import icecube.phys_services
        from I3Tray import I3Tray

        from icecube.simprod.util import BasicCounter
        from icecube.simprod.segments import GenerateCosmicRayMuons, PropagateMuons

        # Instantiate a tray.
        tray = I3Tray()

        summary = dataclasses.I3MapStringDouble()
        summaryfile     = self.GetParameter('summaryfile')
        if os.path.exists(summaryfile): 
           summary = ReadI3Summary(summaryfile)
        tray.context['I3SummaryService'] = summary

        randomService = icecube.phys_services.I3SPRNGRandomService(
            self.seed, self.nproc, self.procnum)\
            if not self.usegslrng else phys_services.I3GSLRandomService(seed = self.seed*self.nproc+self.procnum)

        tray.context["I3RandomService"] = randomService

        tray.AddModule("I3Reader", "read",
                       FilenameList=self.inputfilelist)

        tray.AddSegment(PropagateMuons, "propagator",
                        RandomService=randomService,
                        CylinderLength=self.length,
                        CylinderRadius=self.radius,
                        InputMCTreeName=self.inputmctreename,
                        OutputMCTreeName=self.outputmctreename,
                        PROPOSAL_config_file=self.proposalconfigfile)

        tray.AddModule(BasicCounter, "count_events",
                       Streams=[icecube.icetray.I3Frame.DAQ],
                       name="Generated Events",
                       Stats=stats)


        tray.AddModule("I3Writer", "writer",
                       filename=self.outputfile,
                       Streams=[
                            icecube.icetray.I3Frame.TrayInfo,
                            icecube.icetray.I3Frame.DAQ,
                            icecube.icetray.I3Frame.Stream('S'),
                            icecube.icetray.I3Frame.Stream('M')])

        

        # Execute the tray.
        tray.Execute()
        
        tray.PrintUsage()

        summary = tray.context['I3SummaryService']
        for k in tray.Usage():
            stats[str(k.key())+":usr"] = k.data().usertime
            summary[str(k.key())+":usr"] = k.data().usertime

            stats[str(k.key())+":sys"] = k.data().systime
            summary[str(k.key())+":sys"] = k.data().systime

            stats[str(k.key())+":ncall"] = k.data().ncall
            summary[str(k.key())+":ncall"] = k.data().ncall

        WriteI3Summary(summary, summary_out)


        del tray
        return 0




if __name__ == '__main__':
   import logging
   logging.basicConfig()
   rootLogger = logging.getLogger('')
   rootLogger.setLevel(logging.INFO)

   stats = {}
   prop = MuonPropagator()
   prop.ExecuteOpts(stats)
