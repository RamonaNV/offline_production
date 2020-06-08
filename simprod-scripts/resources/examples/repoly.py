# just load in argparse so the help function can be quickly accessed
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--infile", type=str,
                   help="File in which coincidences will be inserted")
parser.add_argument("-o", "--outfile", type=str,
                    help="Output filename")
parser.add_argument("-s", "--seed", type=int, default=0,
                    help="Seed used in the random service")
parser.add_argument("--backgroundfile", type=str,
                   help="Pre-generated background CORSIKA events")
parser.add_argument("--summaryfile", type=str, default=None,
                    help="(Optional) name of summary file. If not specified, no summary file will be made.")
parser.add_argument("--propagate-background", dest="propagate",
                        default=False, action="store_true", required=False,
                   help="Propagate background events")
parser.add_argument("--mctype", type=str, default="LeptonInjector_numu",
                   help="Primary Type (corsika, nugen, LeptonInjector, etc). To avoid overcounting with corsika.")
parser.add_argument("--MCTreeName", type=str, default="I3MCTree", dest="mctree_name",
                   help="Name of I3MCTree objects to merge")
parser.add_argument("--BGTreeName", type=str, default="BackgroundI3MCTree_preMuonProp", dest="bgtree_name",
                   help="Put coincident events into this separate MCTree. Leave empty to combine.")
parser.add_argument("--TimeWindow", type=float, default=40.0, dest="time_window",
                   help="Coincident event time window in *microseconds*.")
parser.add_argument("--rate", type=float, default=float('nan'),
                   help="Event rate (if NaN polyplopia should get this from background file)")

args = parser.parse_args()


# import time to measure runtime
import time
start_time = time.time()

# do the imports

import icecube
from icecube import dataio, phys_services, polyplopia, dataclasses
from I3Tray import *
from icecube.simprod import segments


# set argparse arguments
infile      = args.infile
outfile     = args.outfile
RNGSeed     = args.seed
SummaryFile = args.summaryfile

MCtype             = args.mctype
MCTreeName         = args.mctree_name
BGTreeName         = args.bgtree_name
TimeWindow         = args.time_window * I3Units.microsecond
Rate               = args.rate



# setup icetray
tray = I3Tray()

randomService = phys_services.I3GSLRandomService(RNGSeed)
tray.context["I3RandomService"] = randomService

tray.context['I3SummaryService'] = dataclasses.I3MapStringDouble()


# Read pre-generated CORSIKA file with CR showers. Assuming the events follow natural spectrum (e.g. Hoerandel)
# Used to inject background CR events from


# Configure tray
tray.AddSegment(dataio.I3Reader, "reader", filenamelist=[infile])

#nugen labels all initial neutrinos before reaching interaction volume as primaries
#other generators can be combined with background to begin with
def check_if_neut_inice(frame,mctree_name="I3MCTree"):
    if "nugen" not in infile.lower():
        return True
    tree = frame[mctree_name]
    most_energy_inice = dataclasses.get_most_energetic_inice(tree)
    if most_energy_inice is not None:
        return True
    else:
        return False
tray.Add(check_if_neut_inice,"check_inice",
         mctree_name = MCTreeName,
         Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

tray.AddSegment(segments.PolyplopiaSegment,"merge",
               mctype = MCtype,
               mctree_name = MCTreeName,
               bgfile = args.backgroundfile,
               separate_coincident_mctree_name = BGTreeName,
               timewindow = TimeWindow,
               rate = Rate*I3Units.hertz,
               RandomService = randomService)

tray.AddModule("Rename","rename_mctree", Keys=[MCTreeName,'SignalI3MCTree'])
tray.AddModule("Rename","rename_mmc", Keys=['MMCTrackList','SignalMMCTrackList'])

if args.propagate:
    tray.AddSegment(segments.PropagateMuons, 'propagator',
               RandomService= randomService,
               SaveState=True,
               InputMCTreeName=BGTreeName,
               OutputMCTreeName="BackgroundI3MCTree",
#               **PROPOSALParams
) 

def CombineSignal(frame,signal_name='SignalI3MCTree',bg_name='BackgroundI3MCTree'):

    signal_tree = frame[signal_name]
    combined_tree = dataclasses.I3MCTree(frame[bg_name])
    polyplopia.MergeMCTrees(combined_tree, signal_tree,0)
    frame["I3MCTree"] = combined_tree

    signal_mmc = None
    if 'MMCTrackList' in frame and 'SignalMMCTrackList' in frame:
        combined_mmc = frame['MMCTrackList']
        signal_mmc = frame['SignalMMCTrackList']
        polyplopia.MergeMMCInfo(combined_mmc, signal_mmc,0)
        del frame['MMCTrackList']
        frame['MMCTrackList'] = combined_mmc
    elif 'SignalMMCTrackList' in frame:
        signal_mmc = frame['SignalMMCTrackList']
        frame['MMCTrackList'] = signal_mmc


tray.Add(CombineSignal,"combine_trees",
             Streams=[icetray.I3Frame.DAQ])


#tray.AddModule("Rename","rename_mmcbg", Keys=['MMCTrackList','ignalMMCTrackList'])



tray.Add("I3Writer", Filename=outfile)

# run the tray
tray.Execute()
tray.Finish()

if SummaryFile:
    from icecube.simprod import util
    
    summary = tray.context['I3SummaryService']
    util.WriteI3Summary(summary, SummaryFile)
    

end_time = time.time()
print("done")
print("That took "+str(end_time - start_time)+" seconds.")
