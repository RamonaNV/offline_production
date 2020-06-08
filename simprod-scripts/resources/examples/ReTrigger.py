#!/usr/bin/env python

from icecube import icetray, dataclasses
from os.path import expandvars

@icetray.traysegment
def ReTrigger(tray, name,
    RunID = None,
    GCDFile = None,
    DeleteKeys = [],
    TimeShiftSkipKeys=[],
    FilterTrigger=True):

    """
    Read photon-propagated (MCPE) files, simulate noise, PTM response, DOMLaunches, and trigger.

    :param RandomService: the name of a random service to be used by the tank response
    :param RunID: Number of run that will be writtend to I3EventHeader
    :param GCDFile: path to GCD file to read first
    :param TimeShiftSkipKeys: list of keys that should be time-shifted. Default: shift all Time-like objects.
    :param FilterTrigger: remove events that don't pass any trigger.
    """
    
    from icecube import icetray, dataclasses, dataio, phys_services
    from icecube import trigger_sim
    from icecube import simclasses,recclasses
    from icecube import lilliput, common_variables

    from I3Tray import I3Units


    tray.Add("Delete", Keys=["I3TriggerHierarchy"]+DeleteKeys)
    tray.Add("Rename", Keys=["TimeShift","TimeShiftOrig"])

        
    timeshiftargs={'SkipKeys':TimeShiftSkipKeys}
    tray.AddSegment(trigger_sim.TriggerSim,
                    name+'_triggersim',
                    gcd_file=dataio.I3File(GCDFile), # for trigger auto-configuration
                    run_id = RunID,
                    prune = True,
                    time_shift = True,
                    time_shift_args = timeshiftargs,
                    filter_mode = FilterTrigger
                    )

    def TotalTimeShift(f):
        if "TimeShiftOrig" in f:
            tshift = f["TimeShift"].value+f["TimeShiftOrig"].value
            del f["TimeShift"]
            del f["TimeShiftOrig"]
            f["TimeShift"] = dataclasses.I3Double(tshift)
        return True
    tray.Add(TotalTimeShift,Streams=[icetray.I3Frame.DAQ] )
    

if __name__ == "__main__":
	from optparse import OptionParser
	from I3Tray import *
	from icecube import dataio

	usage = "usage: %prog [options]"
	parser = OptionParser(usage)

	parser.add_option("-o", "--outfile",
                  default="test_flashes.i3", 
                  dest="OUTFILE", 
                  help="Write output to OUTFILE (.i3{.gz} format)")
	parser.add_option("-i", "--infile",
                  default="", 
                  dest="INFILE", 
                  help="Read input from INFILE (.i3{.gz} format)")
	parser.add_option("-r", "--runid",
                  type=int,
                  dest="RUNID", 
                  help="Run Id")
	parser.add_option("-g", "--gcd",
                  default=expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55380_corrected.i3.gz"),
                  dest="GCDFILE", 
                  help="Read geometry from GCDFILE (.i3{.gz} format)")

	(options,args) = parser.parse_args()
	tray = I3Tray()
	tray.AddModule("I3Reader", "reader",filenamelist=[options.GCDFILE]+[options.INFILE])

	def is_Q(frame):
	    return frame.Stop==frame.DAQ

	simulation_keeps = [
            "BadDomsList",
            "BadDomsListSLC",
            "I3Calibration",
            "I3DetectorStatus",
            "I3Geometry",
            'BackgroundI3MCTreePEcounts',
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
            'I3MCWeightDict',
            'LeptonInjectorProperties',
            'MCPrimary',
            'MCPrimaryInfo',
            'MMCTrackList',
            'PolyplopiaInfo',
            'PolyplopiaPrimary',
            'RNGState',
            'WIMP_params', # Wimp-sim
            'InIceRawData', # keep raw data for now
            'IceTopRawData',
           ]


	tray.AddModule("Keep", "keep_before_retrigger",
                   keys = simulation_keeps,
                   If=is_Q
                   )


	tray.AddSegment(ReTrigger,"retrigger", 
		GCDFile = options.GCDFILE,
		TimeShiftSkipKeys=["I3DST11Header","I3DST12Header","I3DST11","I3DST12",
		               "I3DST13Header","I3DST16Header","I3DST13","I3DST16",
		               "OfflineIceTopHLCPulseInfo",
		               "PoleMuonLlhFitFitParams",
		               "PoleL2BayesianFitFitParams",
		               "PoleL2IpdfGConvolute_2itFitParams",
		               "PoleL2MPEFitFitParams",
		               "I3Calibration",
		               "I3DetectorStatus",
		               "I3Geometry",
		               "FilterMask",
		               "BadDomsListSLC",
		               "BadDomsList",
		               "TankPulseMergerExcludedStations",
		               "CorsikaWeightMap",
		               "ClusterCleaningExcludedStations",
		               "OfflineIceTopHLCPulseInfo",
		               "PoleMuonLlhFitFitParams",
		               ],
		RunID = options.RUNID,
		)
	tray.AddModule("I3Writer", Filename = options.OUTFILE,
            Streams=[icetray.I3Frame.Simulation, icetray.I3Frame.DAQ, icetray.I3Frame.TrayInfo] )
	tray.Execute()
	tray.Finish()
