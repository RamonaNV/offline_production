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
print(" =====      Polyplopia      =====")
print(" ===== Insert coincidences  =====")


# just load in argparse so the help function can be quickly accessed
from argparse import ArgumentParser
from utils import str2bool

parser = ArgumentParser()
parser.add_argument("-i", "--infile",
                    type=str, required=True,
                    help="File in which coincidences will be inserted")

parser.add_argument("-o", "--outfile",
                    type=str, required=True,
                    help="Output filename")

parser.add_argument("-s", "--seed", type=int, required=True,
                    help="Seed used in the random service")

parser.add_argument("--backgroundfile",
                    type=str, required=True,
                    help="Pre-generated background CORSIKA events")

parser.add_argument("--summaryfile",
                    type=str, default=None,
                    help="(Optional) name of summary file. If not specified, no summary file will be made.")

parser.add_argument("--is-corsika", dest="is_corsika",
                    type=str2bool,  default=False,
                    nargs='?', const=True,
                    help="Input file is corsika simulation. Set to avoid overcounting with corsika.",)

parser.add_argument("--MCTreeName", dest="mctree_name",
                    type=str, required=True,
                    help="Name of I3MCTree objects to merge")

parser.add_argument("--OutputMCTreeName", dest="output_mctree_name",
                    type=str, required=True,
                    help="Name of I3MCTree object after Polyplopia")

parser.add_argument("--TimeWindow", dest="time_window",
                    type=float, required=True,
                    help="Coincident event time window in *microseconds*.")

parser.add_argument("--log-level", dest="log_level",
                    type=str, default="WARN",
                    help="Sets the icetray logging level (ERROR, WARN, INFO, DEBUG, TRACE)")

args = parser.parse_args()
print("Called with:")
for key, v in vars(args).items():
    print("{:20s}: {}".format(key, v))
print("")

# import time to measure runtime
import time
start_time = time.time()

print("Importing files... ", end="")
# icecube imports
import icecube
from icecube import dataio, phys_services, polyplopia, dataclasses, icetray
from I3Tray import *
print("done")

# set argparse arguments
infile      = args.infile
outfile     = args.outfile
RNGSeed     = args.seed
bgfile      = args.backgroundfile
SummaryFile = args.summaryfile

if args.is_corsika:
    MCtype = "corsika"
else:
    MCtype = "not_corsika"

MCTreeName  = args.mctree_name
TimeWindow  = args.time_window * I3Units.microsecond
Rate        = float("nan") # get from background file
log_level   = args.log_level

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

randomService = phys_services.I3GSLRandomService(RNGSeed)
tray.context["I3RandomService"] = randomService

tray.context['I3SummaryService'] = dataclasses.I3MapStringDouble()

# Read pre-generated CORSIKA file with CR showers. Assuming the events follow natural spectrum (e.g. Hoerandel)
# Used to inject background CR events from
background = polyplopia.CoincidentI3ReaderService(bgfile)

# Configure tray
tray.AddSegment(dataio.I3Reader, "reader", filenamelist=[infile])

# Ensure that there is at least 1 particle marked as InIce in the I3MCTree:
# NuGen sometimes injects natural rate interactions with no particle marked InIce.
# For these events polyplopia will fail to merge CR coincidences.
def atleast1_inice(frame):
    tree = frame[MCTreeName]
    is_inice = lambda p: p.location_type == dataclasses.I3Particle.InIce
    all_inice = tree.get_filter(is_inice)
    return len(all_inice) > 0

tray.AddModule("PoissonMerger", "merge",
               PrimaryType = MCtype,
               CoincidentEventService = background,
               MCTreeName = MCTreeName,
               SeparateMCTree = "",
               TimeWindow = TimeWindow,
               Rate = Rate,
               RandomServiceName = "I3RandomService",
               If=atleast1_inice)

if args.mctree_name != args.output_mctree_name:
    tray.AddModule("Rename","bookkeeping_mctreenames",
                   Keys=[args.mctree_name, args.output_mctree_name])

tray.Add("I3Writer", Filename=outfile,
            DropOrphanStreams=[icetray.I3Frame.TrayInfo],
            Streams=[icetray.I3Frame.TrayInfo, 
                     icetray.I3Frame.Simulation, 
                     icetray.I3Frame.DAQ, 
                     icetray.I3Frame.Physics])

print("done")
print("Executing Tray...", end="")
tray.Execute()
tray.Finish()
print("done")

if SummaryFile:
    from icecube.simprod import util
    summary = tray.context['I3SummaryService']
    util.WriteI3Summary(summary, SummaryFile)

end_time = time.time()
print("That took "+str(end_time - start_time)+" seconds.")
