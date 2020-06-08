#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/icetray-start
#METAPROJECT /data/user/bsmithers/metaprojects/snobo/py3-v4.0.1/

# Ben Smithers
# benjamin.smithers@mavs.uta.edu

# a simple cascade injector using the Simple Injector
# For creating the generation-level events for the Chiba collaboration meeting

# August 2019

from optparse import OptionParser

#Create parser to accept arguments from user
parser = OptionParser()
parser.add_option("-s","--seed",dest="seed",default=1,type="int",   
                  help="Require a unique number under 4e9 for each job.")
parser.add_option("-e","--energy",dest="energy",default="100000",type=float,
                  help="Energy at which to inject particles")
parser.add_option("-o","--outfile",dest="outfile",type="string",    
                  help="Outfile name",
                  default = None)
parser.add_option("-n","--nEvents",dest="nEvents",default=100,type="int",
                  help="Number of frames to generate")

options,args = parser.parse_args()

print("Importing IceCube modules")
# import icecube stuff
from I3Tray import *
from icecube import dataclasses, phys_services, dataio
from icecube import simple_injector, sim_services
from icecube.icetray import I3Units

seed         = int(options.seed)
try:
    energy   = float(options.energy)*I3Units.GeV
except ValueError:
    raise Exception("Can't cast the energy, {}, as a float".format(energy))
nEvents      = int(options.nEvents)
outfile      = options.outfile

print("Building the IceTray")
# construct a simple tray that 1. creates the DAQ frames, and 2. fills them with an MC tree
tray = I3Tray()
tray.AddService("I3GSLRandomServiceFactory", Seed = seed)
tray.AddModule("I3InfiniteSource",Stream = icetray.I3Frame.DAQ)
tray.AddModule("I3SimpleInjector","injector", 
                ParticleType    = 11, #pdg code for an electron (want cascades)
                EnergyMin       = energy,
                EnergyMax       = energy,
                # particles will be put inside a cylinder of these parameters
                ZMin            = -100.,
                ZMax            = 100.,
                InjectionRadius = 500.,
                # location of the cylinder, [x,y,z], in icecube coordinates
                # note IceCube extends +/- 500 meters in Z
                CylCenter       = [0,0,-400*I3Units.meter],
                GammaIndex      = 2) # not really used since the particles are monoenergetic 

tray.AddModule("I3Writer",Filename=outfile)
print("Executing the IceTray")
tray.Execute(nEvents)
