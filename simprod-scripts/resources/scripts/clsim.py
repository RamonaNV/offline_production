#!/usr/bin/env python

"""
 GPU Photon propagation
"""

import os
from os.path import expandvars
from icecube.simprod.util import CombineHits, DrivingTime, SetGPUEnvironmentVariables
from icecube.simprod.util import ReadI3Summary, WriteI3Summary
from icecube.simprod.util.fileutils import download,untar,isurl
from icecube.simprod.util import simprodtray, arguments
from icecube.simprod.util.simprodtray import RunI3Tray
import argparse
from icecube import clsim
from icecube import polyplopia
from icecube.production_histograms import ProductionHistogramModule
from icecube.production_histograms.histogram_modules.simulation.mcpe_module import I3MCPEModule
from icecube.simprod import segments
from icecube import icetray, dataclasses, dataio, phys_services, interfaces
from I3Tray import I3Tray, I3Units


def add_args(parser):
    """
    Args:
        parser (argparse.ArgumentParser): the command-line parser
    """
    arguments.add_gcdfile(parser)
    arguments.add_inputfilelist(parser)
    arguments.add_outputfile(parser)
    arguments.add_summaryfile(parser)
    arguments.add_enablehistogram(parser)
    arguments.add_histogramfilename(parser)

    arguments.add_nproc(parser)
    arguments.add_procnum(parser)
    arguments.add_seed(parser)  
    arguments.add_usegslrng(parser)  

    arguments.add_icemodellocation(parser)
    arguments.add_icemodel(parser)
    arguments.add_holeiceparametrization(parser)
    arguments.add_oversize(parser)
    arguments.add_efficiency(parser)

    arguments.add_proposalparams(parser)
    arguments.add_propagatemuons(parser, True)

    arguments.add_photonseriesname(parser)
    
    arguments.add_gpu(parser)
    arguments.add_usegpus(parser, False)

    parser.add_argument("--no-RunMPHitFilter", dest="runmphitfilter",
                        default=True, action="store_false", required=False,
                        help="Don't run polyplopia's mphitfilter")
    parser.add_argument("--MCTreeName", dest="mctreename",
                        default="I3MCTree", type=str, required=False,
                        help="Name of MCTree frame object")
    parser.add_argument("--RawPhotonSeriesName", dest="rawphotonseriesname",
                        default=None, type=str, required=False,
                        help="Raw Photon Series Name")
    parser.add_argument("--UseGeant4", dest="usegeant4",
                        default=False, action="store_true", required=False,
                        help="Enable Geant4 propagation")
    parser.add_argument("--no-KeepMCTree", dest="keepmctree",
                        default=True, action="store_false", required=False,
                        help='Delete propagated MCTree')
    parser.add_argument("--keep-pure-background", dest="keeppurebackground",
                        default=False, action="store_false", required=False,
                        help='Keep events without primary hits')


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
    if params['gpu'] is not None and params['usegpus']:
        SetGPUEnvironmentVariables(params['gpu'])

    if params['propagatemuons']:
        if params['usegslrng']:
           randomServiceForPropagators = phys_services.I3GSLRandomService(seed=2 * (params['seed'] * params['nproc'] + params['procnum']))
        else:
           randomServiceForPropagators = phys_services.I3SPRNGRandomService(seed=params['seed'], nstreams=params['nproc'] * 2, streamnum=params['nproc'] + params['procnum'])
        tray.context['I3PropagatorRandomService'] = randomServiceForPropagators

        tray.AddModule("Rename", "rename_corsika_mctree",
                       Keys=[params['mctreename'], params['mctreename'] + '_preMuonProp'])
        tray.AddSegment(segments.PropagateMuons, 'propagator',
                        RandomService=randomServiceForPropagators,
                        **params['proposalparams'])

    tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
                    GCDFile=params['gcdfile'],
                    RandomService=tray.context['I3RandomService'],
                    UseGPUs=params['usegpus'],
                    UseCPUs=not params['usegpus'],
                    IceModelLocation=os.path.join(params['icemodellocation'], params['icemodel']),
                    UnshadowedFraction=params['efficiency'],
                    UseGeant4=params['usegeant4'],
                    DOMOversizeFactor=params['oversize'],
                    MCTreeName=params['mctreename'],
                    MCPESeriesName=params['photonseriesname'],
                    PhotonSeriesName=params['rawphotonseriesname'],
                    HoleIceParameterization=params['holeiceparametrization'])

    if params['runmphitfilter']:
        tray.AddModule("MPHitFilter", "hitfilter",
                       HitOMThreshold=1,
                       RemoveBackgroundOnly=not params['keeppurebackground'],
                       I3MCPESeriesMapName=params['photonseriesname'])

    if params['enablehistogram'] and params['histogramfilename']:
        tray.AddModule(ProductionHistogramModule,
                       Histograms=[I3MCPEModule],
                       OutputFilename=params['histogramfilename'])

    if not params['keepmctree']:
        logger.info("discarding %s" % (params['photonseriesname']))
        tray.Add("Delete", "clean_mctruth",
                 Keys=[params['mctreename'], params['mctreename'] + '_preSampling'])


def main():
    """
     GPU Photon propagation
    """
    # Get Params
    parser = argparse.ArgumentParser(description="ClSim script")
    add_args(parser)
    params = vars(parser.parse_args())  # dict()

    # Process Params
    inputfilenamelist = [params['gcdfile']] + params['inputfilelist']
    
    # Execute Tray
    summary = RunI3Tray(params, configure_tray, "ClSim",
                        summaryfile=params['summaryfile'],
                        inputfilenamelist=inputfilenamelist,
                        outputfile=params['outputfile'],
                        seed=params['seed'],
                        nstreams=params['nproc'],
                        streamnum=params['procnum'],
                        usegslrng=params['usegslrng'])


if __name__ == "__main__":
    main()
