#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT: /data/user/nwandkowsky/tarballs/simulation.V05-01-02/build/simulation.V05-01-02

# import required icecube-related stuff
from icecube import icetray, dataclasses, simclasses, dataio
from icecube.icetray import I3Units
from I3Tray import I3Tray

# command line options required to configure the simulation
from optparse import OptionParser
from os.path import expandvars

#./photon_syst.py -s 1 -r 1 -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_IC86_Merged.i3.gz -o /data/user/nwandkowsky/data/test/NuE.i3.zst --flavor NuE -n 10 --gamma 1.5 --energy_min 5 --energy_max 100000 --usegpu False


usage = "usage: %prog [options] "
#outputfile"
parser = OptionParser(usage)
parser.add_option("-g", "--gcd",default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_IC86_Merged.i3.gz",
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345, dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1, dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=100, dest="NUMEVENTS", help="The number of events per run")
parser.add_option("--flavor", choices=('NuMu', 'NuE', 'NuTau'), default='NuE', dest="FLAVOR", help="Flavor of neutrinos to generate")
parser.add_option("--gamma", type="float", default=2., dest="GAMMA", help="Power law index to use for generation")
parser.add_option("--energy_min", type="float", default=1., dest="EMIN", help="minimal energy in units of TeV")
parser.add_option("--energy_max", type="float", default=10., dest="EMAX", help="maximal energy in units of TeV")      
parser.add_option("-o","--output", default="generator_output.i3", dest="OUTPUT", help="output file name")
parser.add_option("--icemodel", default="Spice3", dest="ICEMODEL", help="icemodel name")
parser.add_option("-d","--domos", type="float", default="5.", dest="DOMOS", help="dom oversizing")
parser.add_option("--domeff", type="float", default="1.17", dest="DOMEFF", help="dom efficiency")
parser.add_option("--usecpu", action="store_false", default="False", dest="USECPU", help="use cpu?")
parser.add_option("--holeice",default="as_50", dest="HOLEICE", help="Holeice file")
parser.add_option("-i","--input", default="generator_output.i3", dest="INPUT", help="output file name")
# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

outfile = options.OUTPUT

print("Command line options parsed...")

import os, sys

tray = I3Tray()

# import phys_services which includes rng
from icecube import phys_services

print("RNG being initiated...")
# set up a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 100000000,
    streamnum = options.RUNNUMBER)
    
# and one for the propagators, that should be hit slightly less often
prandomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 200000000,
    streamnum = options.RUNNUMBER+100000000)

tray.context['I3RandomService'] = randomService

tray.AddModule('I3Reader', 'reader',  filenamelist = [options.GCDFILE, options.INPUT])

tray.Add('Delete','initial_clean',Keys=['PhotonSeriesMap'])

# 3. propagate secondary particles and determine losses
from icecube.simprod.segments.PropagateMuons import PropagateMuons
tray.AddSegment(PropagateMuons,"PropagateMuons", RandomService=prandomService)

# ice model options
if options.ICEMODEL == "SpiceMie":
        clsimIceModel = expandvars("$I3_SRC/clsim/resources/ice/spice_mie")
elif options.ICEMODEL == "SpiceLea":
        clsimIceModel = expandvars("$I3_SRC/clsim/resources/ice/spice_lea")
elif options.ICEMODEL == "spice3_2":
        clsimIceModel = expandvars("$I3_SRC/ice-models/resources/models/spice_3.2")
elif options.ICEMODEL == "spice3_2_scat":
        clsimIceModel = expandvars("$I3_SRC/ice-models/resources/models/spice_3.2_scat")
elif options.ICEMODEL == "spice3_2_abs":
        clsimIceModel = expandvars("$I3_SRC/ice-models/resources/models/spice_3.2_abs")
elif options.ICEMODEL == "spice3_2_scatabs":
        clsimIceModel = expandvars("$I3_SRC/ice-models/resources/models/spice_3.2_scatabs")
else:
        raise RuntimeError("Unknown ice model: %s", IceModel)

# 4. propagate photons
from icecube import clsim
tray.Add(clsim.I3CLSimMakePhotons,
                       UseCPUs=False,
                       UseGPUs=True,
                       MCTreeName="I3MCTree",
                       OutputMCTreeName=None,
                       MMCTrackListName="MMCTrackList",
                       PhotonSeriesName="PhotonSeriesMap",
                       ParallelEvents=100,
                       TotalEnergyToProcess=100.*I3Units.TeV,
                       RandomService=tray.context['I3RandomService'],
                       IceModelLocation=clsimIceModel,
                       UseCascadeExtension=False,
                       StopDetectedPhotons=True,
                       DoNotParallelize=False,
                       DOMOversizeFactor=options.DOMOS,
                       UnshadowedFraction=options.DOMEFF,
                       HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.nominal"),
                       WavelengthAcceptance=None,
                       )

holeiceparameterization=expandvars("$I3_SRC/ice-models/resources/models/holeice_msu/"+options.HOLEICE)

from icecube import clsim
tray.AddSegment(clsim.I3CLSimMakeHitsFromPhotons,"makePhotons",
            PhotonSeriesName="PhotonSeriesMap",
            MCPESeriesName="I3MCPESeriesMap",
            RandomService='I3RandomService',
            DOMOversizeFactor=options.DOMOS,
            UnshadowedFraction=(0.99),
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

# clean up unnecessary keys
tray.Add('Delete', Keys=['PhotonSeriesMap', 'I3MCTree', 'I3MCTree_sliced', 'MMCTrackList'])

# dump content of frame to see if relevant frame objects are there
#tray.AddModule('Dump', "dumpy")

# 5. write some output data
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
