#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT /data/user/nwandkowsky/tarballs/icerec.V05-01-07

# ./L2.py -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz -i det_muongun.i3.zst -o L2_muongun.i3.zst
# ./L2.py -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz -i det_nugen.i3.zst -o L2_nugen.i3.zst

from I3Tray import *
from icecube import dataio, phys_services
#from icecube import icetray, dataclasses

import time

start_time = time.asctime()
print 'Started:', start_time

# handling of command line arguments  
from optparse import OptionParser
parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)

parser.add_option("-i", "--input", action="store", type="string", default="", dest="INPUT", help="Output i3 file")
parser.add_option("-o", "--output", action="store", type="string", default="", dest="OUTPUT", help="Output i3 file")
parser.add_option("-g", "--gcd", action="store", type="string", default="", dest="GCD", help="GCD file for input i3 file")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

outfile = options.OUTPUT

tray = I3Tray()

tray.AddModule('I3Reader', 'reader', FilenameList=[options.GCD, options.INPUT])
print "Reading input file...", options.INPUT
        
if 'I3_DATA' in os.environ:
    TableDir = os.path.expandvars("$I3_DATA/photon-tables/")
elif os.path.isdir('/cvmfs/icecube.opensciencegrid.org/data'):
    TableDir = os.path.expandvars("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/")
else:
    raise Exception('cannot find I3_DATA or cvmfs, no photon tables available')
SplineDir = os.path.join(TableDir, "splines")
SplineRecoAmplitudeTable = os.path.join(SplineDir, 'InfBareMu_mie_abs_z20a10_V2.fits')
SplineRecoTimingTable = os.path.join(SplineDir, 'InfBareMu_mie_prob_z20a10_V2.fits')

from icecube.filterscripts.all_filters import OnlineFilter

tray.AddSegment(OnlineFilter, "OnlineFilter",
                    sdstarchive=False, # our input data is SDST archive data
                    SplineRecoAmplitudeTable = SplineRecoAmplitudeTable,
                    SplineRecoTimingTable = SplineRecoTimingTable,                    
                    simulation=True, decode = False, If=lambda f: True,
                    slop_split_enabled = False,
                    vemcal_enabled=False, 
                    gfu_enabled=False,
                    needs_wavedeform_spe_corr = False,
                    hese_followup=False,
                    estres_followup=False,
                    ehealert_followup=False
                    )

# make random service
seed = os.getpid()
filter_mask_randoms = phys_services.I3GSLRandomService(seed)

# override MinBias Prescale
from icecube.filterscripts import filter_globals
filterconfigs = filter_globals.filter_pairs + filter_globals.sdst_pairs
print(filterconfigs)

# Generate filter Masks for all P frames
from icecube import filter_tools
tray.AddModule(filter_tools.FilterMaskMaker, "MakeFilterMasks",
               OutputMaskName = filter_globals.filter_mask,
               FilterConfigs = filterconfigs,
               RandomService = filter_mask_randoms)

# Merge the FilterMasks
tray.AddModule("OrPframeFilterMasks", "make_q_filtermask",
               InputName = filter_globals.filter_mask,
               OutputName = filter_globals.qfilter_mask)

#Q+P frame specific keep module needs to go first, as KeepFromSubstram
#will rename things, let's rename post keep.  
def is_Q(frame):
    return frame.Stop==frame.DAQ

gcd_keeps = ['I3Geometry','I3Calibration','I3DetectorStatus','BadDomsList','BadDomsListSLC','BadDomsListHLC']
simulation_keeps = [
        'BackgroundI3MCTree',
        'BackgroundI3MCTreePEcounts',
        'BackgroundI3MCPESeriesMap',
        'BackgroundI3MCTree_preMuonProp',
        'BackgroundMMCTrackList',
        'BeaconLaunches',
        'CorsikaInteractionHeight',
        'CorsikaWeightMap',
        'EventProperties',
        'GenerationSpec',
        'I3LinearizedMCTree',
        'I3MCTree',
        'I3MCTreePEcounts',
        'I3MCTree_preMuonProp',
        'I3MCPESeriesMap',
        'I3MCPulseSeriesMap',
        'I3MCPulseSeriesMapParticleIDMap',
        'I3MCWeightDict',
        'LeptonInjectorProperties',
        'MCHitSeriesMap',
        'MCPrimary',
        'MCTrack',
        'MCPrimaryInfo',
        'MMCTrackList',
	'PolyplopiaCount',
        'PolyplopiaInfo',
        'PolyplopiaPrimary',
        'RNGState',
        'SignalI3MCPEs',
        'SimTrimmer', # for SimTrimmer flag
        'TimeShift', # the time shift amount
        'WIMP_params', # Wimp-sim
       ]

keep_before_merge = filter_globals.q_frame_keeps + [
                        'InIceDSTPulses', # keep DST pulse masks
                        'IceTopDSTPulses',
                        'CalibratedWaveformRange', # keep calibration info
                        'UncleanedInIcePulsesTimeRange',
                        'SplitUncleanedInIcePulses',
                        'SplitUncleanedInIcePulsesTimeRange',
                        'SplitUncleanedInIceDSTPulsesTimeRange',
                        'CalibrationErrata',
                        'SaturationWindows',
                        'InIceRawData', # keep raw data for now
                        'IceTopRawData',
                       ] + simulation_keeps + gcd_keeps

tray.AddModule("Keep", "keep_before_merge",
               keys = keep_before_merge,
               If=is_Q
               )


## second set of prekeeps, conditional on filter content, based on newly created Qfiltermask
#Determine if we should apply harsh keep for events that failed to pass any filter
##  Note: excluding the sdst_streams entries

tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheckAll",
               FilterNameList = filter_globals.filter_streams,
               FilterResultName = filter_globals.qfilter_mask,
               DecisionName = "PassedAnyFilter",
               DiscardEvents = False, 
               Streams = [icetray.I3Frame.DAQ]
               )
def do_save_just_superdst(frame):
    if frame.Has("PassedAnyFilter"):
        if not frame["PassedAnyFilter"].value:
            return True    #  <- Event failed to pass any filter.  
        else:
            return False # <- Event passed some filter

    else:
        print("Failed to find key frame Bool!!")
        return False

keep_only_superdsts = filter_globals.keep_nofilterpass+[
                         'PassedAnyFilter',
                         'InIceDSTPulses',
                         'IceTopDSTPulses',
                         'SplitUncleanedInIcePulses',
                         'SplitUncleanedInIcePulsesTimeRange',
                         'SplitUncleanedInIceDSTPulsesTimeRange',
                         'RNGState',
                         ] + simulation_keeps + gcd_keeps
tray.AddModule("Keep", "KeepOnlySuperDSTs",
               keys = keep_only_superdsts,
               If = do_save_just_superdst
               )

## Now clean up the events that not even the SuperDST filters passed on.
tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheckSDST",
               FilterNameList = filter_globals.sdst_streams,
               FilterResultName = filter_globals.qfilter_mask,
               DecisionName = "PassedKeepSuperDSTOnly",
               DiscardEvents = False,
               Streams = [icetray.I3Frame.DAQ]
               )

def dont_save_superdst(frame):
    if frame.Has("PassedKeepSuperDSTOnly") and frame.Has("PassedAnyFilter"):
        if frame["PassedAnyFilter"].value:
            return False  #  <- these passed a regular filter, keeper
        elif not frame["PassedKeepSuperDSTOnly"].value:
            return True    #  <- Event failed to pass SDST filter.  
        else:
            return False # <- Event passed some  SDST filter
    else:
        print("Failed to find key frame Bool!!")
        return False

tray.AddModule("Keep", "KeepOnlyDSTs",
               keys = filter_globals.keep_dst_only
                      + ["PassedAnyFilter","PassedKeepSuperDSTOnly",
                         filter_globals.eventheader] + gcd_keeps,
               If = dont_save_superdst
               )
   
## Frames should now contain only what is needed.  now flatten, write/send to server
## Squish P frames back to single Q frame, one for each split:
tray.AddModule("KeepFromSubstream","null_stream",
               StreamName = filter_globals.NullSplitter,
               KeepKeys = filter_globals.null_split_keeps,
               )

from icecube.phys_services.which_split import which_split
in_ice_keeps = filter_globals.inice_split_keeps + filter_globals.onlinel2filter_keeps
in_ice_keeps = in_ice_keeps + ['I3EventHeader',
                               'SplitUncleanedInIcePulses',
                               'SplitUncleanedInIcePulsesTimeRange',
                               'TriggerSplitterLaunchWindow',
                               'I3TriggerHierarchy',
                               'GCFilter_GCFilterMJD'] + gcd_keeps
tray.AddModule("Keep", "inice_keeps",
               keys = in_ice_keeps,
               If = which_split(split_name=filter_globals.InIceSplitter),
               )


tray.AddModule("KeepFromSubstream","icetop_split_stream",
               StreamName = filter_globals.IceTopSplitter,
               KeepKeys = filter_globals.icetop_split_keeps,
               )

# Apply small keep list (SuperDST/SmallTrig/DST/FilterMask for non-filter passers
# Remove I3DAQData object for events not passing one of the 'filters_keeping_allraw'
tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheck",
               FilterNameList = filter_globals.filters_keeping_allraw,
               FilterResultName = filter_globals.qfilter_mask,
               DecisionName = "PassedConventional",
               DiscardEvents = False,
               Streams = [icetray.I3Frame.DAQ]
               )

## Clean out the Raw Data when not passing conventional filter
def I3RawDataCleaner(frame):
    if not (('PassedConventional' in frame and 
             frame['PassedConventional'].value == True) or 
            ('SimTrimmer' in frame and
             frame['SimTrimmer'].value == True)
           ):
        frame.Delete('InIceRawData')
        frame.Delete('IceTopRawData')

tray.AddModule(I3RawDataCleaner,"CleanErrataForConventional",
               Streams=[icetray.I3Frame.DAQ])
    
# DO L2 PROCESSING HERE    

from icecube.filterscripts.offlineL2.level2_all_filters import OfflineFilter
tray.AddSegment(OfflineFilter, "OfflineFilter",
        dstfile=False,
        mc=True,
        doNotQify=True,
        photonicsdir=TableDir
        )

def streamcut(frame):
    return frame["I3EventHeader"].sub_event_stream!='IceTopSplit' and frame["I3EventHeader"].sub_event_stream!='NullSplit'
tray.Add(streamcut)

tray.Add('Delete', 'final_delete', Keys=[
  'ClusterCleaningExcludedTanks', 'IceTopDSTPulses', 'IceTopPulses', 
  'CleanIceTopRawData', 'IceTopRawData',
  'OfflineIceTopHLCTankPulses', 'OfflineIceTopHLCVEMPulses', 'OfflineIceTopSLCVEMPulses',
  'SimTrimmer',
  'TankPulseMergerExcludedTanks', 
  'I3MCPESeriesMap', 'I3MCPESeriesMap_081','I3MCPESeriesMap_090','I3MCPESeriesMap_095','I3MCPESeriesMap_099','I3MCPESeriesMap_108','I3MCPESeriesMap_117',
  #'InIceRawData'
  ])

#tray.Add('Dump')

tray.AddModule('I3Writer', 'writer',
    DropOrphanStreams=[icetray.I3Frame.DAQ],
    Streams=[  icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    filename=outfile)

tray.AddModule('TrashCan', 'thecan')

tray.Execute()
tray.Finish()

del tray

stop_time = time.asctime()

print 'Started:', start_time
print 'Ended:', stop_time
