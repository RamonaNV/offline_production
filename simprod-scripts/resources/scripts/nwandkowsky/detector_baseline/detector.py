#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT: /data/user/nwandkowsky/tarballs/simulation.V05-01-02/build/simulation.V05-01-02

# import required icecube-related stuff
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units
from I3Tray import I3Tray

# command line options required to configure the simulation
from optparse import OptionParser
from os.path import expandvars

#./mcpe_nugen.py -s 2 -r 1 -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_IC86_Merged.i3.gz --domos 5 --domeff 0.99 --holeice flasher_p1=0.3_p2=0.0 -i /data/ana/Cscd/StartingEvents/NuGen_new/NuE/medium_energy/photon_spice3_2/1/photon_00000001.i3.zst -o nue_me_mcpe.i3.bz2 -b /data/ana/Cscd/StartingEvents/CORSIKA_bg/12531/photon_spice3_2/1/photon_00000001.i3.zst


usage = "usage: %prog [options] "
#outputfile"
parser = OptionParser(usage)
parser.add_option("-g", "--gcd",default="/home/nwandkowsky/workspace/data/GCD/GeoCalibDetectorStatus_IC86.55697_V2.i3",
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-i", "--input", action="store", type="string", default="", dest="INPUT", help="Input file to process")
parser.add_option("-o", "--output", action="store", type="string", default="", dest="OUTPUT", help="Output i3 file")    
parser.add_option("-s", "--seed",type="int",default=12345, dest="SEED", help="Initial seed for the random number generator")                  
parser.add_option("-r", "--runnumber", type="int", default=1, dest="RUNNUMBER", help="The run number for this simulation")   
parser.add_option("--holeice",default="as_50", dest="HOLEICE", help="Holeice file")
parser.add_option("--domos", type="float", default=1., dest="DOMOS", help="dom oversizing parameter")
parser.add_option("--domeff", type="float", default=1., dest="DOMEFF", help="dom efficiency parameter")
parser.add_option("-b", "--bgfile", action="store", type="string", default="", dest="BGFILE", help="Output i3 file") 
# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

GCD = options.GCDFILE
infile = options.INPUT
outfile = options.OUTPUT

print("Command line options parsed...")

holeiceparameterization=expandvars("$I3_SRC/ice-models/resources/models/holeice_msu/"+options.HOLEICE)
               
import os, sys

tray = I3Tray()
tray.AddModule("I3Reader", "reader", filenamelist=[GCD, infile] )

# import phys_services which includes rng
from icecube import phys_services
# set up a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 100000000,
    streamnum = options.RUNNUMBER)
tray.context['I3RandomService'] = randomService

prandomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 200000000,
    streamnum = options.RUNNUMBER+100000000)

tray.Add('Delete',Keys=["I3MCTree","MMCTrackList"])

from icecube.simprod.segments.PropagateMuons import PropagateMuons
tray.AddSegment(PropagateMuons,"PropagateMuons", RandomService=prandomService)

def RemoveInvisible(frame):
    return len(frame["PhotonSeriesMap"])>0;
tray.Add(RemoveInvisible)(("Streams",[icetray.I3Frame.DAQ]));

from icecube import polyplopia
tray.AddService("CoincidentI3ReaderServiceFactory","BackgroundService")(("FileName",options.BGFILE));

tray.AddModule("PoissonPEMerger")(
    ("BaseIsBackground",False),
    ("CoincidentEventService","BackgroundService"),
    ("MCTreeName","I3MCTree_preMuonProp"),
    ("RandomService","I3RandomService"),
    ("PhotonsToMerge","PhotonSeriesMap"),
    ("MMCTrackName","MMCTrackList")
    );

tray.AddModule("I3GeometryDecomposer", "decomposeGeometry")

from icecube import clsim
tray.AddSegment(clsim.I3CLSimMakeHitsFromPhotons,"makePhotons",
            PhotonSeriesName="PhotonSeriesMap",
            MCPESeriesName="I3MCPESeriesMap", 
            RandomService='I3RandomService', 
            DOMOversizeFactor=options.DOMOS,
            UnshadowedFraction=(1.17),
            HoleIceParameterization=holeiceparameterization
        )

# Delete all MCPEs we're not operating on
def empty_mcpe(frame):
    entries = 0
    for k in frame.keys():
        if isinstance(frame[k], simclasses.I3MCPESeriesMap):
            entries = entries + len(frame[k])
    return entries>0
tray.AddModule(empty_mcpe, Streams=[icetray.I3Frame.DAQ])
        
map_in = 'I3MCPESeriesMap'
map_out = 'I3MCPESeriesMap_'+str(options.DOMEFF)
domeff = float(options.DOMEFF)/1.17
print map_in, map_out, domeff

# Downsampling from input pe series (map_in with highest dom eff = 1.17) to desired output pe series (map_out)
tray.AddModule("I3DownsampleMCPE", InputName=map_in, OutputName=map_out, SampleFrac = domeff, RandomService = 'I3RandomService')

from icecube.simprod import segments
tray.AddSegment(segments.DetectorSim, "DetectorSim",
    RandomService = 'I3RandomService',
    RunID = options.RUNNUMBER,
    GCDFile = GCD,
    KeepMCHits = False,
    KeepPropagatedMCTree = True,
    KeepMCPulses = False,
    SkipNoiseGenerator = False,
    LowMem = True,
    InputPESeriesMapName = map_out
    )

# acquire some more MC truth information
def GetNeutrino(frame):
    def sanitize(particle):
        if particle is None:
            return dataclasses.I3Particle()
        else:
            return particle

    mcTree = frame["I3MCTree"]
    primary = None
    neutrino = None
    for p in mcTree:
        if mcTree.depth(p) != 0: continue

        if p.is_neutrino:
            if neutrino is None or p.energy > neutrino.energy:
                neutrino = p

        if primary is None or p.energy > primary.energy:
            primary = p
    del frame["MCPrimary"]
    frame["MCPrimary"] = sanitize(primary)
tray.AddModule(GetNeutrino, "GetNeutrino", Streams=[icetray.I3Frame.DAQ])


def GetMCTrack(frame):
    
   # Get the track from an MCTree. If there is no track, get the hadronic cascade.
    mcTree = frame["I3MCTree"]

    trackParticle = None
    cascadeParticle = None
    numCascades = 0
    neutrino = None

    for p in mcTree:
        depth = mcTree.depth(p)
        if depth == 0:
            # if not p.is_neutrino:
            #     raise RuntimeError("primary particle is not a neutrino!")

            if neutrino is None or p.energy > neutrino.energy:
                neutrino = p

        if depth != 1: continue # depth==0 is the root (assume it is the primary neutrino)

        if p.type in [dataclasses.I3Particle.ParticleType.MuPlus,
                      dataclasses.I3Particle.ParticleType.MuMinus,
                      dataclasses.I3Particle.ParticleType.TauPlus,
                      dataclasses.I3Particle.ParticleType.TauMinus]:
           # if trackParticle is not None:
           #     raise RuntimeError("got multiple leptons from a single neutrino.")

            trackParticle = p
        else:
            if cascadeParticle is None or p.energy > cascadeParticle.energy:
                cascadeParticle = p
            numCascades += 1

    theTrack = None

    if trackParticle is not None:
        theTrack = trackParticle
    else:
        if numCascades == 0: theTrack = None
        if numCascades == 1: theTrack = cascadeParticle
        if neutrino is None:
            raise RuntimeError("Internal error. Cascades found, but no neutrino in MCTree.")
        theTrack = neutrino

    if theTrack is None:
        raise RuntimeError("no MC track could be found in MCTree")

    # shift the vertex to the point of closest approach to the origin (0,0,0)
    a = - (theTrack.pos.x*theTrack.dir.x + theTrack.pos.y*theTrack.dir.y + theTrack.pos.z*theTrack.dir.z)
    newPos = dataclasses.I3Position(theTrack.pos.x + theTrack.dir.x * a,
                                    theTrack.pos.y + theTrack.dir.y * a,
                                    theTrack.pos.z + theTrack.dir.z * a)
    newTime = theTrack.time + a/dataclasses.I3Constants.c

    # generate a "reconstructed" particle from the MCTrack
    outputTrack = dataclasses.I3Particle()
    outputTrack.shape = dataclasses.I3Particle.ParticleShape.InfiniteTrack
    outputTrack.pos = newPos
    outputTrack.dir = theTrack.dir
    outputTrack.time = newTime
    outputTrack.fit_status = dataclasses.I3Particle.FitStatus.OK
    outputTrack.location_type = dataclasses.I3Particle.LocationType.InIce

    frame["MCTrack"] = outputTrack
    
tray.AddModule(GetMCTrack, "GetMCTrack", Streams=[icetray.I3Frame.DAQ])

#tray.AddModule('Dump', "dumpy")

tray.Add('Delete',Keys=[map_out,'I3MCPESeriesMap','I3MCPulseSeriesMapPrimaryIDMap','IceTopRawData','PhotonSeriesMap'])


tray.AddModule('I3Writer', 'writer',
    Streams=[icetray.I3Frame.Stream('S'), icetray.I3Frame.TrayInfo, icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    filename=outfile)

# final clean-up and execution
tray.AddModule('TrashCan', 'YesWeCan')
print("Executing...")
tray.Execute()
print("Finish!")
tray.Finish()

del tray
