from .. import ipmodule
from ..util import ReadI3Summary, WriteI3Summary
from ..util import CombineHits, DrivingTime, choose_max_efficiency
import os.path


class MuonGunGenerator(ipmodule.ParsingModule):
    """IceProd module for ``MuonGun`` simulations

    It defines parameters for and executes
    :func:`icecube.simprod.segments.GenerateCosmicRayMuons` and
    :func:`icecube.simprod.segments.PropagateMuons`.

    **Module parameters:**

        :attr:`nproc`
            Number of processes for RNG
        :attr:`pronum`
            Process number
        :attr:`seed`
            RNG seed
        :attr:`gcdfile`
            GeoCalibDetStatus filename
        :attr:`nevents`
            Number of generated events
        :attr:`model`
            Primary cosmic-ray flux parametrization
        :attr:`gamma`
            Power law spectral index
        :attr:`offset`
            Power law offset in GeV
        :attr:`emin`
            Mininum generated energy in GeV
        :attr:`emax`
            Maximum generated energy in GeV
        :attr:`length`
            Cylinder length in m
        :attr:`radius`
            Cylinder radius in m
        :attr:`x`
            Cylinder x-position in m
        :attr:`y`
            Cylinder y-position in m
        :attr:`z`
            Cylinder z-position in m
        :attr:`length_dc`
            Inner cylinder length in m
        :attr:`radius_dc`
            Inner cylinder radius in m
        :attr:`x_dc`
            Inner cylinder x-position in m
        :attr:`y_dc`
            Inner cylinder y-position in m
        :attr:`z_dc`
            Inner cylinder z-position in m
        :attr:`deepcore`
            Use inner cylinder
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

        self.AddParameter("gcdfile", "GeoCalibDetStatus filename", "")
        self.AddParameter("nevents", "number of generated events", 1)
        self.AddParameter("model", "primary cosmic-ray flux parametrization",
                          "Hoerandel5_atmod12_SIBYLL")
        self.AddParameter("gamma", "power law spectral index", 2.)
        self.AddParameter("offset", "power law offset in GeV", 700.)
        self.AddParameter("emin", "mininum generated energy in GeV", 1e4)
        self.AddParameter("emax", "maximum generated energy in GeV", 1e7)
        self.AddParameter("length", "cylinder length in m", 1600.)
        self.AddParameter("radius", "cylinder radius in m", 800.)
        self.AddParameter("x", "cylinder x-position in m", 0.)
        self.AddParameter("y", "cylinder y-position in m", 0.)
        self.AddParameter("z", "cylinder z-position in m", 0.)
        self.AddParameter("length_dc", "inner cylinder length in m", 500.)
        self.AddParameter("radius_dc", "inner cylinder radius in m", 150.)
        self.AddParameter("x_dc", "inner cylinder x-position in m", 46.3)
        self.AddParameter("y_dc", "inner cylinder y-position in m", -34.9)
        self.AddParameter("z_dc", "inner cylinder z-position in m", -300.)
        self.AddParameter("deepcore", "use inner cylinder", False)
        self.AddParameter('propagate_muons','Run PROPOSAL.', True)
        self.AddParameter('propagate_photons','Run ClSim.', True)

        self.AddParameter("outputfile", "output filename", "corsika.i3")
        self.AddParameter("summaryfile", "Summary filename", "muongun.json")

        self.AddParameter('HistogramFilename', 'Histogram filename.', None)
        self.AddParameter('EnableHistogram', 'Write a SanityChecker histogram file.', False)
        self.AddParameter("natural_rate", "Sample natural rate muon bundles", False)
        self.AddParameter("GPU", 
                          "Graphics Processing Unit number (shoud default to environment if None)",
                          None)
        self.AddParameter("UseGPUs", "Use Graphics Processing Unit",False)
        self.AddParameter("UseOnlyDeviceNumber", "Use only this device.", 0)
        self.AddParameter("oversize","over-R: DOM radius oversize scaling factor",5)
        self.AddParameter("holeiceparametrization", "Location of hole ice param files", 
                          os.path.expandvars("$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm"))
        self.AddParameter("efficiency","overall DOM efficiency correction",[1.00])
        self.AddParameter("volumecyl","set volume to regular cylinder (set to False for 300m spacing from the DOMs)",True)
        self.AddParameter("IceModelLocation","Location of ice model param files", os.path.expandvars("$I3_BUILD/ice-models/resources/models")) 
        self.AddParameter("IceModel","ice model subdirectory", "spice_3.2") 
        self.AddParameter("PhotonSeriesName","Photon Series Name","I3MCPESeriesMap") 
        self.AddParameter("RawPhotonSeriesName","Raw Photon Series Name",None) 
        self.AddParameter('KeepMCTree','Delete propagated MCTree otherwise',True)
        self.AddParameter("UseGSLRNG","Use I3GSLRandomService",False) 


    def Execute(self, stats):
        if not ipmodule.ParsingModule.Execute(self, stats):
            return 0

        import icecube.icetray
        import icecube.dataclasses
        import icecube.dataio
        import icecube.phys_services
        from I3Tray import I3Tray

        from ..util import BasicCounter
        from ..segments import GenerateCosmicRayMuons, PropagateMuons, GenerateNaturalRateMuons



        if self.gpu is not None and self.usegpus:
           os.putenv("CUDA_VISIBLE_DEVICES",str(self.gpu))
           os.putenv("COMPUTE",":0."+str(self.gpu))
           os.putenv("GPU_DEVICE_ORDINAL",str(self.gpu))


        # Instantiate a tray.
        tray = I3Tray()

        summary = icecube.dataclasses.I3MapStringDouble()
        tray.context['I3SummaryService'] = summary

        randomService = icecube.phys_services.I3SPRNGRandomService(
            self.seed, self.nproc, self.procnum)\
            if not self.usegslrng else icecube.phys_services.I3GSLRandomService(self.seed)
        tray.context["I3RandomService"] = randomService

        tray.AddModule("I3InfiniteSource","TheSource",
                       Prefix=self.gcdfile,
                       Stream=icecube.icetray.I3Frame.DAQ)

        if self.natural_rate:
           tray.AddSegment(GenerateNaturalRateMuons, "muongun",
                        NumEvents=self.nevents,
                        mctree_name="I3MCTree_preMuonProp",
                        flux_model="GaisserH4a_atmod12_SIBYLL",
                        )


        else:
           # Configure tray segment that actually does stuff.
           tray.AddSegment(GenerateCosmicRayMuons, "muongun",
                        mctree_name="I3MCTree_preMuonProp",
                        num_events=self.nevents,
                        flux_model=self.model,
                        gamma_index=self.gamma,
                        energy_offset=self.offset,
                        energy_min=self.emin,
                        energy_max=self.emax,
                        cylinder_length=self.length,
                        cylinder_radius=self.radius,
                        cylinder_x=self.x,
                        cylinder_y=self.y,
                        cylinder_z=self.z,
                        inner_cylinder_length=self.length_dc,
                        inner_cylinder_radius=self.radius_dc,
                        inner_cylinder_x=self.x_dc,
                        inner_cylinder_y=self.y_dc,
                        inner_cylinder_z=self.z_dc,
                        use_inner_cylinder=self.deepcore)

        if self.propagate_muons:
            tray.AddSegment(PropagateMuons, "propagator",
                            RandomService=randomService,
                            CylinderLength=self.length,
                            CylinderRadius=self.radius,
                            InputMCTreeName="I3MCTree_preMuonProp",
                            OutputMCTreeName="I3MCTree")

        tray.AddModule(BasicCounter, "count_events",
                       Streams=[icecube.icetray.I3Frame.DAQ],
                       name="Generated Events",
                       Stats=stats)

        if self.propagate_photons:
            if not self.propagate_muons: 
                raise BaseException("You have to propagate muons if you want to propagate photons")

            if type(self.efficiency) == list or type(self.efficiency) == tuple:
                if len(self.efficiency) == 1:
                    efficiency=float(self.efficiency[0])
                elif len(self.efficiency) > 1:
                    efficiency=map(float,self.efficiency)
                elif len(self.efficiency) > 1:
                    raise Exception("Configured empty efficiency list")
            else:
                efficiency = choose_max_efficiency(self.efficiency)

            from icecube import clsim
            try:                
                tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
                                GCDFile = self.gcdfile,
                                RandomService = randomService,
                                UseGPUs = self.usegpus,
                                UseOnlyDeviceNumber = self.useonlydevicenumber,
                                UseCPUs= not self.usegpus,
                                IceModelLocation = os.path.join(self.icemodellocation,self.icemodel),
                                UnshadowedFraction = efficiency,
                                UseGeant4 = False,
                                DOMOversizeFactor = self.oversize,
                                MCTreeName = "I3MCTree", 
                                MCPESeriesName = self.photonseriesname,
                                PhotonSeriesName = self.rawphotonseriesname,
                                HoleIceParameterization = self.holeiceparametrization)

                from icecube import polyplopia
                tray.AddModule("MPHitFilter","hitfilter",
                               HitOMThreshold=1,
                               RemoveBackgroundOnly=False,
                               I3MCPESeriesMapName=self.photonseriesname)
                
            except AttributeError as e:
                print(e)
                print("Nevermind...not propagating photons.")


        if not self.keepmctree:
            self.logger.info("discarding %s" % (self.photonseriesname))
            tray.Add("Delete","clean_mctruth",
                     Keys=["I3MCTree",'I3MCTree_preSampling'])


        if self.enablehistogram and self.histogramfilename:         
            from icecube.production_histograms import ProductionHistogramModule
            from icecube.production_histograms.histogram_modules.simulation.mctree_primary import I3MCTreePrimaryModule
            from icecube.production_histograms.histogram_modules.simulation.mctree import I3MCTreeModule
            from icecube.production_histograms.histogram_modules.simulation.mcpe_module import I3MCPEModule
        
            tray.AddModule(ProductionHistogramModule, 
                           Histograms = [I3MCTreePrimaryModule, I3MCTreeModule,I3MCPEModule],
                           OutputFilename = self.histogramfilename)

        tray.AddModule("I3Writer", "writer",
                       filename=self.outputfile,
                       streams=[
                            icecube.icetray.I3Frame.TrayInfo,
                            icecube.icetray.I3Frame.DAQ, 
                            icecube.icetray.I3Frame.Stream('S'), 
                            icecube.icetray.I3Frame.Stream('M')]
                       )


        # Execute the tray.
        tray.Execute()
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

        WriteI3Summary(summary, self.summaryfile)

        del tray
        return 0
