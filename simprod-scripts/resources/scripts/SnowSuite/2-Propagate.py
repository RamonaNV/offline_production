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
print(" ===== Propagate Particles  =====")

# just load in argparse so the help function can be quickly accessed
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-s", "--seed", dest="seed",
                    default=1, type=int,
                    help="Require a unique number under 4e9 for each job.")

parser.add_argument("-i", "--infile", dest="infile",
                    default=None, type=str, required=True, 
                    help="the full path to the input file")

parser.add_argument("-o","--outfile",dest="outfile",
                    default=None, type=str, required=True,
                    help="the full path to the output file")

parser.add_argument("-g", "--gcdfile", dest="gcdfile",
                    default=None, type=str, required=True,
                    help="the full path to the calibration file")

parser.add_argument("--log-level", dest="log_level",
                    type=str, default="WARN",
                    help="Sets the icetray logging level (ERROR, WARN, INFO, DEBUG, TRACE)")

args = parser.parse_args()

print("")
print("Called with:")
for key, v in vars(args).items():
    print("{:10s}: {}".format(key, v))
print("")

# import time to measure runtime
import time
start_time = time.time()

# pass the arguments to python variables
seed        = int(args.seed)
infile      = args.infile
outfile     = args.outfile
gcdfile     = args.gcdfile
log_level   = args.log_level


print("Importing necessary packages...", end="")
# import python packages
import os
from os.path import expandvars
import numpy as np

# import icecube ones
from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube import phys_services
from icecube.icetray import I3Units
from icecube.simprod import segments
print("done")

# set icetray logging level
log_levels = {"error" : icetray.I3LogLevel.LOG_ERROR,
              "warn" : icetray.I3LogLevel.LOG_WARN,
              "info" : icetray.I3LogLevel.LOG_INFO,
              "debug" : icetray.I3LogLevel.LOG_DEBUG,
              "trace" : icetray.I3LogLevel.LOG_TRACE}
if log_level.lower() in log_levels.keys():
    icetray.set_log_level(log_levels[log_level.lower()])


print("Preparing Tray...", end="")
tray = I3Tray()

randomServicePropagators = phys_services.I3GSLRandomService(seed=seed)

tray.Add(dataio.I3Reader, "reader",
         filenamelist=[gcdfile, infile])

tray.AddModule("Rename", keys=["I3MCTree", "I3MCTree_preMuonProp"])

tray.AddSegment(segments.PropagateMuons, "PropagateMuons",
                RandomService = randomServicePropagators)

tray.AddModule("Delete", keys=["I3MCTree_preMuonProp"])

tray.AddModule("I3Writer", Filename=outfile,
               DropOrphanStreams=[icetray.I3Frame.TrayInfo],
               Streams=[icetray.I3Frame.TrayInfo,
                        icetray.I3Frame.Simulation,
                        icetray.I3Frame.DAQ,
                        icetray.I3Frame.Physics])

print("done")
print("Executing Tray...", end="")
tray.Execute()
tray.Finish()

end_time = time.time()
print("done")
print("That took "+str(end_time - start_time)+" seconds.")
