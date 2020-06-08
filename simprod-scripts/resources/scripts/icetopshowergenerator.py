#!/usr/bin/env python

"""
Wrapper runs GenerateIceTopShowers Tray segment
"""

import os
from I3Tray import I3Tray, I3Units
from icecube.simprod import segments
from icecube.simprod.util import simprodtray, arguments
from icecube.simprod.util import ReadI3Summary, WriteI3Summary
from icecube import icetray, dataclasses, dataio
from icecube.simprod.util.simprodtray import RunI3Tray
import argparse
import icecube.icetray
import math
from icecube import phys_services, sim_services
from icecube.simprod.util import PrintContext


def add_args(parser):
    """
    Args:
        parser (argparse.ArgumentParser): the command-line parser
    """
    arguments.add_gcdfile(parser)
    arguments.add_inputfilelist(parser)
    arguments.add_outputfile(parser)
    arguments.add_summaryfile(parser)

    arguments.add_nproc(parser)
    arguments.add_procnum(parser)
    arguments.add_seed(parser)
    arguments.add_usegslrng(parser)

    arguments.add_propagatemuons(parser, True)
    
    parser.add_argument("--samples", dest="samples",
                        default=1, type=int, required=False,
                        help='Number of samples of the same CORSIKA shower')
    parser.add_argument("--RunId", dest="runid",
                        default=None, type=str, required=False,
                        help='RunId (or dataset)')
    parser.add_argument("--x", dest="x",
                        default=0., type=float, required=False,
                        help='Core is place in a disk around (x,y)')
    parser.add_argument("--y", dest="y",
                        default=0., type=float, required=False,
                        help='Core is place in a disk around (x,y)')
    parser.add_argument("--r", dest="r",
                        default=0., type=float, required=False,
                        help='The radius of the disk whithin wich the core is randomly chosen (usually energy-dependent, something like 800 + 600*(log10(E/GeV) - 5) meters')
    parser.add_argument("--unthin-r", dest="unthin_r",
                        default=0.*I3Units.meter, type=float, required=False,
                        help="Radius of sampling region for CORSIKA unthinning. Usually energy-dependent, something like -693.4 + 360.4*log10(E/GeV) -60.4*(log10(E/GeV))^2 +3.34*(log10(E/GeV))^3 (hey, don't ask me...)")                    
    parser.add_argument("--dump", dest="dump",
                        default=False, action="store_true", required=False,
                        help="Dump frame contentes to screen")
    parser.add_argument("--tank-response", dest="tank_response",
                        default='g4', type=str, required=False,
                        help='g4 or param. Use Geant4 or parametrized propagation')
    parser.add_argument("--tank-sampling-radius", dest="tank_sampling_radius",
                        default=25.*I3Units.meter, type=float, required=False,
                        help='default 25 m (for standard MC, for thinned, should be larger or equal than UnThinRadius)')
    parser.add_argument("--raise-observation-level", dest="raise_observation_level",
                        default=0., type=float, required=False,
                        help="Tweak the altitude (in meters) where corsika particles are injected (just in case the corsika observation plane is below the top of the snow)")


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
    # first check some variables that might need to be set
    dataset = 0
    if params['runid'] is None:
       runid = int(10**math.ceil(math.log10(params['nproc']))*dataset + params['procnum'])
    else:
       runid = params['runid']

    # THE THING
    tray.AddSegment(segments.GenerateIceTopShowers, "GenerateIceTopShowers",
                    NSamples = params['samples'],
                    Files = params['inputfilelist'],
                    GCDFile = params['gcdfile'],
                    x=params['x'], y=params['y'], r=params['r']*I3Units.meter, 
                    TankResponse=params['tank_response'],
                    TankSamplingRadius = params['tank_sampling_radius'],
                    UnthinRadius=params['unthin_r'],
                    RunID=runid,
                    RaiseObservationLevel=params['raise_observation_level']*I3Units.m)

    if params['propagatemuons']:
       # segments.PropagateMuons requires a random service, I would prefer to use the same one used above.
       # I could also set the state of randomService to the state of tray.context['I3RandomService'] before running this segment.
       if params['usegslrng']:
           randomService = phys_services.I3GSLRandomService(seed=params['seed']*params['nproc']+params['procnum'])
       else:
           randomService = phys_services.I3SPRNGRandomService(seed=params['seed'], nstreams=params['nproc'], streamnum=params['procnum'])

       tray.Add('Rename', Keys=['I3MCTree', 'I3MCTree_preMuonProp'])
       tray.AddSegment(segments.PropagateMuons, "PropagateMuons",
                       RandomService = randomService)
    
    tray.AddModule(PrintContext,"ctx")
    if params['dump']:
       tray.AddModule("Dump","dump")
    
    
def main():
    """
    Wrapper runs GenerateIceTopShowers Tray segment
    """
    # Get Params
    parser = argparse.ArgumentParser(description="IceTopShowerGenerator script")
    add_args(parser)
    params = vars(parser.parse_args()) # dict()
    
    # Execute Tray
    summary = RunI3Tray(params, configure_tray, "IceTopShowerGenerator",
                        summaryfile=params['summaryfile'], 
                        summaryin=dataclasses.I3MapStringDouble(),
                        outputfile=params['outputfile'],
                        outputstreams=[icetray.I3Frame.DAQ],
                        seed=params['seed'],
                        nstreams=params['nproc'],
                        streamnum=params['procnum'],
                        usegslrng=params['usegslrng'])


if __name__ == "__main__":
    main()
