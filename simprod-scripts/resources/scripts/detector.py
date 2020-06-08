#!/usr/bin/env python
"""
Runs the detector simulation for IC86 for primary events.

Optionally simulates coincident background in-ice using PolyPlopia.
The number of coincidences may be random or forced, in the latter
case the event is weighted according to the chance probability for
such a multiplicity.

Minimum inputs:
 - gcdfile: A GCD file
 - inputfile: an input file with I3MCPEs

Optional inputs:
 - summaryfile: A XML summary file
 - BackGroundMCFile: a file with background events.

By default, background simulation is off. To turn it on, either set
  EnablePolyplopia = True
or
  NumberOfPrimaries = <integer larger than 1>

If PolyPlopia is used, you can define the PolyplopiaRate = <some rate>,
or set it to zero, then the rate is computed from summaryfile.
"""

import os
from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, simclasses, dataio
from icecube.icetray import I3Frame
from icecube.simprod import segments
from icecube.simprod.util import ReadI3Summary, WriteI3Summary
from icecube.simprod.util import simprodtray, arguments
from icecube.simprod.util.simprodtray import RunI3Tray
import argparse
from icecube import clsim
from icecube import phys_services
from icecube import sim_services
from icecube import vuvuzela
from icecube import DOMLauncher
from icecube import trigger_sim
from icecube.production_histograms import ProductionHistogramModule
from icecube.production_histograms.histogram_modules.simulation.pmt_response import PMTResponseModule
from icecube.production_histograms.histogram_modules.simulation.dom_mainboard_response import InIceResponseModule
from icecube.production_histograms.histogram_modules.simulation.trigger import TriggerModule
from icecube.production_histograms.histograms.simulation.noise_occupancy import NoiseOccupancy


def add_args(parser):
    """
    Args:
        parser (argparse.ArgumentParser): the command-line parser
    """
    arguments.add_gcdfile(parser)
    arguments.add_inputfile(parser)
    arguments.add_outputfile(parser)
    arguments.add_summaryfile(parser)
    arguments.add_enablehistogram(parser)
    arguments.add_histogramfilename(parser)

    arguments.add_nproc(parser)
    arguments.add_procnum(parser)
    arguments.add_seed(parser)
    arguments.add_usegslrng(parser)

    parser.add_argument("--MCType", dest="mctype",
                        default='corsika_weighted', type=str, required=False,
                        help='Generator particle type')
    parser.add_argument("--UseLinearTree", dest="uselineartree",
                        default=False, action="store_true", required=False,
                        help='Use I3LinearizedMCTree for serialization')
    parser.add_argument("--MCPrescale", dest="mcprescale",
                        default=100, type=int, required=False,
                        help='Prescale for keeping additional Monte Carlo info in the frame')
    parser.add_argument("--IceTop", dest="icetop",
                        default=False, action="store_true", required=False,
                        help='Do IceTop Simulation?')
    parser.add_argument("--Genie", dest="genie",
                        default=False, action="store_true", required=False,
                        help='Assume separate Genie MCPEs and BG MCPEs')
    parser.add_argument("--no-FilterTrigger", dest="filtertrigger",
                        default=True, action="store_false", required=False,
                        help="Don't filter untriggered events")
    parser.add_argument("--no-Trigger", dest="trigger",
                        default=True, action="store_false", required=False,
                        help="Don't run trigger simulation")
    parser.add_argument("--LowMem", dest="lowmem",
                        default=False, action="store_true", required=False,
                        help='Low Memory mode')
    parser.add_argument("--no-BeaconLaunches", dest="beaconlaunches",
                        default=True, action="store_false", required=False,
                        help="Don't simulate beacon launches")
    parser.add_argument("--TimeShiftSkipKeys", dest="timeshiftskipkeys",
                        default=[], type=arguments.str_comma_list, required=False,
                        help='Skip keys in the triggersim TimeShifter')
    parser.add_argument("--SampleEfficiency", dest="sampleefficiency",
                        default=0.0, type=float, required=False,
                        help='Resample I3MCPESeriesMap for different efficiency')
    parser.add_argument("--GeneratedEfficiency", dest="generatedefficiency",
                        default=0.0, type=float, required=False,
                        help='Generated efficiency for resampling')
    parser.add_argument("--RunID", dest="runid",
                        default=0, type=int, required=False,
                        help='Run ID')
    parser.add_argument("--MCPESeriesName", dest="mcpeseriesname",
                        default='I3MCPESeriesMap', type=str, required=False,
                        help='Name of MCPESeriesMap in frame')
    parser.add_argument("--DetectorName", dest="detectorname",
                        default='IC86', type=str, required=False,
                        help='Name of detector')
    parser.add_argument("--SkipKeys", dest="skipkeys",
                        default=[], type=arguments.str_comma_list, required=False,
                        help='Skip keys for the writer')


def configure_tray(tray, params, stats, logger):
    """
    Configures the I3Tray instance: adds modules, segments, services, etc.

    Args:
        tray (I3Tray): the IceProd tray instance
        params (dict): command-line arguments (and default values)
                            referenced as dict entries; see add_args()
        stats (dict): dictionary that collects run-time stats
        logger (logging.Logger): the logger for this script
    """
    if params['mcpeseriesname'] != 'I3MCPESeriesMap':
        def mover(fr):
            if 'I3MCPESeriesMap' in fr:
                del fr['I3MCPESeriesMap']
            fr['I3MCPESeriesMap'] = fr[params['mcpeseriesname']]
        tray.Add(mover, "move_MCPESeries", Streams=[icetray.I3Frame.DAQ])

    tray.AddSegment(segments.DetectorSegment, "detector",
                    gcdfile=params['gcdfile'],
                    mctype=params['mctype'],
                    uselineartree=params['uselineartree'],
                    detector_label=params['detectorname'],
                    runtrigger=params['trigger'],
                    filtertrigger=params['filtertrigger'],
                    stats=stats,
                    icetop=params['icetop'],
                    genie=params['genie'],
                    prescale=params['mcprescale'],
                    lowmem=params['lowmem'],
                    BeaconLaunches=params['beaconlaunches'],
                    TimeShiftSkipKeys=params['timeshiftskipkeys'],
                    SampleEfficiency=params['sampleefficiency'],
                    GeneratedEfficiency=params['generatedefficiency'],
                    RunID=params['runid'],
                    KeepMCHits=not params['procnum'] % params['mcprescale'],
                    KeepPropagatedMCTree=not params['procnum'] % params['mcprescale'],
                    KeepMCPulses=not params['procnum'] % params['mcprescale'])

    if params['enablehistogram'] and params['histogramfilename']:
        tray.AddModule(ProductionHistogramModule,
                       Histograms=[PMTResponseModule,
                                   InIceResponseModule,
                                   TriggerModule,
                                   NoiseOccupancy],
                       OutputFilename=params['histogramfilename'])


def main():
    # Get Params
    parser = argparse.ArgumentParser(description="IceCube script")
    add_args(parser)
    params = vars(parser.parse_args())  # dict()

    # Execute Tray
    summary = RunI3Tray(params, configure_tray, "IceCube",
                        summaryfile=params['summaryfile'],
                        inputfilenamelist=[params['gcdfile'], params['inputfile']],
                        outputfile=params['outputfile'],
                        outputskipkeys=params['skipkeys'],
                        seed=params['seed'],
                        nstreams=params['nproc'],
                        streamnum=params['procnum'],
                        usegslrng=params['usegslrng'])


if __name__ == "__main__":
    main()
