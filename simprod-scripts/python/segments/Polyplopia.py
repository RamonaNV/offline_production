from I3Tray import *
from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Frame
from icecube.icetray import traysegment
import os.path 
from os.path import expandvars

def SetMultiplicity(frame,mctreelist=[],weightmap="CorsikaWeightMap"):
    if weightmap  in frame:
        wm = frame[weightmap]
    else:
        wm = dataclasses.I3MapStringDouble()
    multiplicity = 0
    for t in mctreelist:
        multiplicity += len(frame[t].primaries)
    wm["Multiplicity"] = multiplicity
    if weightmap  in frame:
       frame.Delete(weightmap)
    frame.Put(weightmap,wm)
    return True



@traysegment
def PolyplopiaSegment(tray, name,
                    mctype='CORSIKA',
                    RandomService=None,
                    mctree_name = "I3MCTree_preMuonProp",
                    separate_coincident_mctree_name = "", # leave empty to combine
                    bgfile = None,
                    timewindow = 40.*I3Units.microsecond,
                    rate = float('nan'),
                    If=lambda f: True
   ):

      """
         There are three scenarios for polyplopia:
            1. bgfile: We are reading background MC from a separate file and
            injecting events to signal (or weighted background) based on 
            a Poisson distribution within the given time window.

            2. we are generating MuonGun bundles and
            injecting events to signal (or weighted background) based on 
            a Poisson distribution within the given time window.
      """

      from icecube import polyplopia, MuonGun

      #tray.AddModule("ParticleMapper","mapprimary") 
      if bgfile: # merge bg into signal
          background = polyplopia.CoincidentI3ReaderService(bgfile)

      else:
          # Use MuonGun
          # Default: use Hoerandel as a template for generating muons.
          model = MuonGun.load_model("Hoerandel5_atmod12_SIBYLL")
          #model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')


          # Generate bundles (even if not 100 percent accurate).
          model.flux.max_multiplicity = 100

          gamma_index=2.6
          energy_offset=700.
          energy_min=1e4
          energy_max=1e7
          cylinder_length=1600.
          cylinder_radius=800.
          cylinder_x=0.
          cylinder_y=0.
          cylinder_z=0.

          # Default: cylinder aligned with z-axis at detector center as injection
          # surface.
          outsurface_center = dataclasses.I3Position(
               cylinder_x*I3Units.m,
               cylinder_y*I3Units.m,
               cylinder_z*I3Units.m)

          outsurface = MuonGun.Cylinder(
               length=cylinder_length*I3Units.m,
               radius=cylinder_radius*I3Units.m,
               center=outsurface_center)

          generator = MuonGun.NaturalRateInjector(outsurface, model.flux, model.energy)


          background = polyplopia.MuonGunBackgroundService()
          background.set_generator(generator)
          background.set_rng(RandomService)
          background.set_rate(rate)
          background.set_mctree_name(mctree_name)

      tray.AddModule("PoissonMerger","merge",
          CoincidentEventService = background,
          PrimaryType = mctype,
          MCTreeName = mctree_name,
          Rate = rate,
          SeparateMCTree = separate_coincident_mctree_name,
          TimeWindow = timewindow)

      return tray



@traysegment
def PolyplopiaPhotons(tray, name,
                    mctype='CORSIKA',
                    RandomService=None,
                    RandomServiceForPropagators=None,
                    mctree_name = "I3MCTree",
                    bgfile = None,
                    GCDFile = None,
                    IceModel = "SpiceLea",
                    IceModelLocation = os.path.expandvars("$I3_BUILD/clsim/resources/ice"),
                    timewindow = 40.*I3Units.microsecond,
                    rate = float('nan'),
                    GPU = None,
                    UseGPUs = True,
                    DOMOversizeFactor = 5,
                    HoleIceParameterization = expandvars("$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm"),
                    Efficiency = 0.99,
                    PhotonSeriesName = "I3MCPESeriesMap",
                    PROPOSALParams=dict(),
                    UsePPC=False,
                    If=lambda f: True
   ):

        from .. import segments
        from I3Tray import I3Units

        tray.AddSegment(segments.PolyplopiaSegment,"coincify",
               RandomService=RandomService,
               mctype=mctype,
               mctree_name = mctree_name,
               separate_coincident_mctree_name = "BackgroundI3MCTree_preMuonProp",
               bgfile = bgfile,
               timewindow = timewindow,
               rate = rate,
           )

        tray.AddModule("Rename","rename_mmc", Keys=['MMCTrackList','SignalMMCTrackList'])

        if RandomServiceForPropagators is None: 
            RandomServiceForPropagators = RandomService

        tray.AddSegment(segments.PropagateMuons, 'propagator',
               RandomService= RandomServiceForPropagators,
               SaveState=True,
               InputMCTreeName="BackgroundI3MCTree_preMuonProp",
               OutputMCTreeName="BackgroundI3MCTree",
               **PROPOSALParams) 
        
        from icecube.simprod import util
        if GPU is not None:
           util.SetGPUEnvironmentVariables(GPU)

        tray.AddModule("Rename","rename_pes", Keys=[PhotonSeriesName,'SignalI3MCPEs'])

        if UsePPC: # Propagate photons
           tray.AddSegment(segments.PPCTraySegment,"ppc_photons",
                           UseGPUs = UseGPUs,
                           UnshadowedFraction = Efficiency,
                           DOMOversizeFactor = DOMOversizeFactor,
                           IceModelLocation = IceModelLocation.replace('clsim','ppc'),
                           IceModel = IceModel,
                           volumecyl = True,
                           gpulib = 'opencl',
                           InputMCTree="BackgroundI3MCTree",
                           keep_empty_events = True,
                           MCPESeriesName = "BackgroundI3MCPESeriesMap")
        else: # use clsim
            from icecube import clsim
            tray.AddSegment(clsim.I3CLSimMakeHits, "makeBackgroundCLSimHits",
                            GCDFile = GCDFile,
                            RandomService = RandomService,
                            UseGPUs = UseGPUs,
                            UseCPUs= not UseGPUs,
                            IceModelLocation = os.path.join(IceModelLocation,IceModel),
                            UnshadowedFraction = Efficiency,
                            UseGeant4 = False,
                            DOMOversizeFactor = DOMOversizeFactor,
                            MCTreeName = "BackgroundI3MCTree",
                            MCPESeriesName = "BackgroundI3MCPESeriesMap",
                            HoleIceParameterization = HoleIceParameterization)

        from icecube import polyplopia
        if mctype.lower() =='corsika':
           WeightMap="CorsikaWeightMap"
        else:
           WeightMap="I3MCWeightDict"

        tray.AddModule("MPHitFilter","hitfilter",
              HitOMThreshold=1,
              RemoveBackgroundOnly=False,
              I3MCPESeriesMapName="BackgroundI3MCPESeriesMap",
              MCTreeName="BackgroundI3MCTree",
              PruneTree=True,
              Filter=False)

        tray.Add(SetMultiplicity,
             mctreelist=[mctree_name,"BackgroundI3MCTree"],
             weightmap="PolyplopiaInfo",
             Streams=[icetray.I3Frame.DAQ])

        # Add to CorsikaWeightMap because people are used to finding it there
        tray.Add(SetMultiplicity,
             mctreelist=[mctree_name,"BackgroundI3MCTree"],
             weightmap=WeightMap,
             Streams=[icetray.I3Frame.DAQ])

        # Make separate Signal, Background and also a Combined Tree
        if mctree_name != 'SignalI3MCTree':
           tray.AddModule("Rename","rename_mctree", Keys=[mctree_name,'SignalI3MCTree'])
           def merge_trees(f):
              if not 'SignalI3MCTree' in f or not 'BackgroundI3MCTree' in f:
                 raise Exception("Could not find either signal or background tree in frame!")

              import copy
              combined = f['SignalI3MCTree']
              background = copy.deepcopy(f['BackgroundI3MCTree'])
              combined.merge(background)
              f[mctree_name] = combined
           tray.Add(merge_trees, Streams=[icetray.I3Frame.DAQ])

        tray.AddModule("Rename","rename_mmc1", Keys=['MMCTrackList','BackgroundMMCTrackList'])
        tray.AddModule("Rename","rename_mmc2", Keys=['SignalMMCTrackList','MMCTrackList'])

        tray.AddModule("I3CombineMCPE","combinepes",
            InputResponses = ['SignalI3MCPEs','BackgroundI3MCPESeriesMap'],
            OutputResponse = PhotonSeriesName
        )

@traysegment
def PolyplopiaMergePEs(tray, name,
                    mctype='CORSIKA',
                    RandomService=None,
                    mctree_name = "I3MCTree",
                    separate_coincident_mctree_name = "", # leave empty to combine
                    bgfile = None,
                    timewindow = 40.*I3Units.microsecond,
                    rate = float('nan'),
                    If=lambda f: True
   ):

      """
      This segment can be used to merge background that already contains I3MCPEs. We are reading background MC from a separate file and
      injecting events to signal (or weighted background) based on 
      a Poisson distribution within the given time window.

      :param mctype: type of signal that we are injecting backgroun onto. Needed to avoid overcounting background.
      :param RandomService: the name of a random service to be used by the tank response
      :param mctree_name: Name of I3MCTree in background file
      :param separate_coincident_mctree_name: Name of output bg tree. If empty, background will be added to main tree.
      :param bgfile: Name of file containing background with I3MCPEs.
      :param mctype: type of signal that we are injecting backgroun onto. Needed to avoid overcounting background.
      :param timewindow: coincidence time window.
      :param rate: rate of background muons cointained in file.
      """
      if mctype.lower() =='corsika':
           WeightMap="CorsikaWeightMap"
      else:
           WeightMap="I3MCWeightDict"

      tray.AddModule("Rename","rename_pes", Keys=["I3MCPESeriesMap",'SignalI3MCPEs'])

      from icecube.polyplopia import PoissonPEMerger
      tray.Add(PoissonPEMerger,"merge",
            RandomService = RandomService,
            BackgroundFile = bgfile,
            I3MCPESeriesName = "SignalI3MCPEs",
            MCTreeName = "I3MCTree",
            MMCTrackName ="MMCTrackList",
            TimeWindow = timewindow,
            Rate = rate)

      tray.AddModule("MPHitFilter","hitfilter",
              HitOMThreshold=1,
              RemoveBackgroundOnly=False,
              I3MCPESeriesMapName="BackgroundI3MCPESeriesMap",
              MCTreeName="BackgroundI3MCTree",
              PruneTree=True,
              Filter=False)

      tray.Add(SetMultiplicity,
             mctreelist=[mctree_name,"BackgroundI3MCTree"],
             weightmap="PolyplopiaInfo",
             Streams=[icetray.I3Frame.DAQ])

      # Add to CorsikaWeightMap because people are used to finding it there
      tray.Add(SetMultiplicity,
             mctreelist=[mctree_name,"BackgroundI3MCTree"],
             weightmap=WeightMap,
             Streams=[icetray.I3Frame.DAQ])

      tray.AddModule("I3CombineMCPE","combinepes",
            InputResponses = ['SignalI3MCPEs','BackgroundI3MCPESeriesMap'],
            OutputResponse = "I3MCPESeriesMap",
        )

 
