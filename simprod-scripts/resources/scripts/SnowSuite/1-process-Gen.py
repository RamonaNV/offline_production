from __future__ import print_function

print(" ===== Welcome to SnowSuite =====")
print("")
print("                /\ ")
print("           __   \/   __")
print("           \_\_\/\/_/_/ ")
print("             _\_\/_/_ ")
print("            __/_/\_\__ ")
print("           /_/ /\/\ \_\ ")
print("                /\ ")
print("                \/ ")
print("")
print("=== Running Lepton Generation ===")

# just load in argparse so the help function can be quickly accessed
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-s", "--seed", dest="seed",
                    type=int, required=True,
                    help="Require a unique number under 4e9 for each job.")

parser.add_argument("-o", "--outfile", dest="outfile",
                    type=str, required=True,
                    help="Outfile name")

parser.add_argument("--outfile-lic", dest="outfile_lic",
                    type=str, required=True,
                    help="Lic Outfile name")

parser.add_argument("-i", "--index", dest="index",
                    type=float, required=True,
                    help="Spectral index, positive")

parser.add_argument("--Emin", dest="Emin",
                    type=float, required=True,
                    help="Min energy to simulate [GeV]")

parser.add_argument("--Emax", dest="Emax",
                    type=float, required=True,
                    help="Maximum energy to smulate [GeV]")

parser.add_argument("-n", "--nEvents", dest="nEvents",
                    type=int, required=True,
                    help="Number of events")

parser.add_argument("--Zmin", dest="Zmin",
                    type=float, required=True,
                    help="Min zenith in degrees, 90 is horizon")

parser.add_argument("--Zmax", dest="Zmax",
                    type=float, required=True,
                    help="Max zenith in degrees, 180 is core")

parser.add_argument("--radius", dest="radius",
                    type=float, required=True,
                    help="Radius around the origin within which to target events/"\
                         "Radius of the vertical cylinder around the origin within which to place events.")

parser.add_argument("--length", dest="length",
                    type=float, required=True,
                    help="Length of the fixed endcaps add to each end of the distance along which to sample interactions/"\
                         "Height of the vertical cylinder around the origin within which to place events.")

parser.add_argument("-t", "--type", dest="type",
                    type=str, required=True,
                    help="Primary neutrino type (SM only). Electron, Muon, Tau, or All")

parser.add_argument("-c","--interaction", dest="interaction",
                    type=str, required=True,
                    help="What kind of interaction do you want to simulate?")

parser.add_argument("--log-level", dest="log_level",
                    type=str, default="WARN",
                    help="Sets the icetray logging level (ERROR, WARN, INFO, DEBUG, TRACE)")

args = parser.parse_args()

print("")
print("Called with:")
for key, v in vars(args).items():
    print("{:17s}: {}".format(key, v))
print("")

# import time to measure runtime
import time
start_time = time.time()

# pass the arguments to python variables
outfile     = args.outfile
seed        = args.seed
Emin        = args.Emin
Emax        = args.Emax
index       = args.index
nEvents     = args.nEvents
Zmin        = args.Zmin
Zmax        = args.Zmax
radius      = args.radius
length      = args.length
nu_type     = args.type
interaction = args.interaction
log_level   = args.log_level

print("Importing necessary packages...", end="")
# import python packages
import os
from os.path import expandvars
import numpy as np
import sys
import warnings

# import icecube ones
from I3Tray import *
from icecube import icetray, dataio, dataclasses, LeptonInjector
from icecube import phys_services
from icecube.icetray import I3Units
print("done")

# check python version
if sys.version_info[0]>2:
    warnings.warn("Warning! LIC files written in Python3 are currently unreadable!")

# set icetray logging level
log_levels = {"error" : icetray.I3LogLevel.LOG_ERROR,
              "warn" : icetray.I3LogLevel.LOG_WARN,
              "info" : icetray.I3LogLevel.LOG_INFO,
              "debug" : icetray.I3LogLevel.LOG_DEBUG,
              "trace" : icetray.I3LogLevel.LOG_TRACE}
if log_level.lower() in log_levels.keys():
    icetray.set_log_level(log_levels[log_level.lower()])

# takes a string representing a particle, returns a pair of it and its anti-particle.
# also accepts 'all', in which case it returns one of every charged lepton
def get_fs_particle_CC(nu_type):
    if nu_type=='Tau' or nu_type=='tau':
        return([dataclasses.I3Particle.ParticleType.TauMinus, dataclasses.I3Particle.ParticleType.TauPlus])
    elif nu_type=='Muon' or nu_type=='muon':
        return([dataclasses.I3Particle.ParticleType.MuMinus , dataclasses.I3Particle.ParticleType.MuPlus])
    elif nu_type=='electron' or nu_type=='Electron':
        return([dataclasses.I3Particle.ParticleType.EMinus  , dataclasses.I3Particle.ParticleType.EPlus])
    elif nu_type=='all' or nu_type =='All':
        return([dataclasses.I3Particle.ParticleType.TauMinus, dataclasses.I3Particle.ParticleType.TauPlus,
                dataclasses.I3Particle.ParticleType.MuMinus , dataclasses.I3Particle.ParticleType.MuPlus,
                dataclasses.I3Particle.ParticleType.EMinus  , dataclasses.I3Particle.ParticleType.EPlus])
    else:
        raise Exception("I don't recogize what a '{}' is. Try Electron, Muon, Tau, or All.".format(nu_type))

def get_fs_particle_NC(nu_type):
    if nu_type=='Tau' or nu_type=='tau':
        return([dataclasses.I3Particle.ParticleType.NuTau, dataclasses.I3Particle.ParticleType.NuTauBar])
    elif nu_type=='Muon' or nu_type=='muon':
        return([dataclasses.I3Particle.ParticleType.NuMu , dataclasses.I3Particle.ParticleType.NuMuBar ])
    elif nu_type=='electron' or nu_type=='Electron':
        return([dataclasses.I3Particle.ParticleType.NuE  , dataclasses.I3Particle.ParticleType.NuEBar  ])
    elif nu_type=='all' or nu_type =='All':
        return([dataclasses.I3Particle.ParticleType.NuTau, dataclasses.I3Particle.ParticleType.NuTauBar,
                dataclasses.I3Particle.ParticleType.NuMu , dataclasses.I3Particle.ParticleType.NuMuBar,
                dataclasses.I3Particle.ParticleType.NuE  , dataclasses.I3Particle.ParticleType.NuEBar])
    else:
        raise Exception("I don't recogize what a '{}' is. Try Electron, Muon, Tau, or All.".format(nu_type))


print("Setting up injectors...", end="")
# need to prepare a list of injectors with which to instantiate the generator module
# the cross sections and final state particles will all depend on which interaction we're doing.
dis_xs_folder = "/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0"
injector_list = []

if interaction.lower()=='gr':
    raise NotImplementedError("GR interaction not implemented!")
    """
    # in a GR interaciton, you have a NuEBar annihilate with an e-. Need to conserve charge + lepton no.
    party_pairs = [[dataclasses.I3Particle.ParticleType.EMinus  , dataclasses.I3Particle.ParticleType.NuEBar  ],
                   [dataclasses.I3Particle.ParticleType.MuMinus , dataclasses.I3Particle.ParticleType.NuMuBar ],
                   [dataclasses.I3Particle.ParticleType.TauMinus, dataclasses.I3Particle.ParticleType.NuTauBar],
                   [dataclasses.I3Particle.ParticleType.Hadrons , dataclasses.I3Particle.ParticleType.Hadrons ]]
    doubly = os.path.join(dis_xs_folder, "single_xs_gr.fits")
    total  = os.path.join(dis_xs_folder, "total_xs_gr.fits")
    """
elif interaction.lower()=='cc':
    final_states = get_fs_particle_CC(nu_type)
    
    for final_state in final_states:
        is_ranged = (final_state==dataclasses.I3Particle.MuMinus  or
                     final_state==dataclasses.I3Particle.MuPlus   or
                     final_state==dataclasses.I3Particle.TauMinus or
                     final_state==dataclasses.I3Particle.TauPlus)
        
        if int(final_state) > 0:
            doubly = os.path.join(dis_xs_folder, "dsdxdy_nu_CC_iso.fits")
            total  = os.path.join(dis_xs_folder, "sigma_nu_CC_iso.fits")
        elif int(final_state) < 0:
            doubly = os.path.join(dis_xs_folder, "dsdxdy_nubar_CC_iso.fits")
            total  = os.path.join(dis_xs_folder, "sigma_nubar_CC_iso.fits")
        else:
            raise Exception("Don't know which cross section to use for {}".format(final_state))
        
        print("Adding CC injector for {} with ranged mode: {}".format(final_state, is_ranged))
        injector = LeptonInjector.injector(
            NEvents                            = int(nEvents/len(final_states)),
            FinalType1                         = final_state,
            FinalType2                         = dataclasses.I3Particle.ParticleType.Hadrons,
            DoublyDifferentialCrossSectionFile = doubly,
            TotalCrossSectionFile              = total,
            Ranged                             = is_ranged)
        injector_list.append(injector)
elif interaction.lower()=='nc':
    final_states = get_fs_particle_NC(nu_type)
    
    for final_state in final_states:
        if int(final_state) > 0:
            doubly = os.path.join(dis_xs_folder, "dsdxdy_nu_NC_iso.fits")
            total  = os.path.join(dis_xs_folder, "sigma_nu_NC_iso.fits")
        elif int(final_state) < 0:
            doubly = os.path.join(dis_xs_folder, "dsdxdy_nubar_NC_iso.fits")
            total  = os.path.join(dis_xs_folder, "sigma_nubar_NC_iso.fits")
        else:
            raise Exception("Don't know which cross section to use for {}".format(final_state))
        
        print("Adding NC injector for {} with ranged mode: {}".format(final_state, False))
        injector = LeptonInjector.injector(
            NEvents                            = int(nEvents/len(final_states)),
            FinalType1                         = final_state,
            FinalType2                         = dataclasses.I3Particle.ParticleType.Hadrons,
            DoublyDifferentialCrossSectionFile = doubly,
            TotalCrossSectionFile              = total,
            Ranged                             = False)
        injector_list.append(injector)
else:
    raise Exception("I don't know what kind of interaction '{}' is. Try GR, CC, or NC".format(interaction))
print("done")

# initialize icetray
print("Constructing I3Tray...", end="")
tray = I3Tray()

randomService = phys_services.I3GSLRandomService(seed=seed)
tray.context["I3RandomService"] = randomService

tray.AddService("I3EarthModelServiceFactory", "Earth")

tray.AddModule("I3InfiniteSource", "TheSource",
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("MultiLeptonInjector",
               EarthModel      = "Earth",
               Generators      = injector_list,
               MinimumEnergy   = Emin * I3Units.GeV,
               MaximumEnergy   = Emax * I3Units.GeV,
               MinimumZenith   = Zmin * I3Units.deg,
               MaximumZenith   = Zmax * I3Units.deg,
               PowerLawIndex   = index,
               InjectionRadius = radius * I3Units.meter,
               EndcapLength    = length * I3Units.meter,
               CylinderRadius  = radius * I3Units.meter,
               CylinderHeight  = length * I3Units.meter,
               MinimumAzimuth  = 0. * I3Units.deg,
               MaximumAzimuth  = 360. * I3Units.deg,           
               RandomService   = "I3RandomService")

# adds a header to the DAQ frames
event_id = 1
def get_header(frame):
    global event_id
    header                  = dataclasses.I3EventHeader()
    header.event_id         = event_id
    header.run_id           = seed
    event_id                += 1
    frame["I3EventHeader"]  = header

tray.AddModule(get_header, streams=[icetray.I3Frame.DAQ]);

tray.AddModule("InjectionConfigSerializer",
               OutputPath=args.outfile_lic)

tray.AddModule("I3Writer", Filename=outfile,
               Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation,
                        icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

print("done")
print("Executing Tray...")
tray.Execute()
tray.Finish()
print("Done!")

end_time = time.time()
print("That took "+str(end_time - start_time)+" seconds.")
