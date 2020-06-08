#!/usr/bin/env python

"""
 GPU Photon propagation
"""
import argparse

from I3Tray import I3Tray, I3Units
from icecube.simprod.util import ReadI3Summary, WriteI3Summary
from icecube.simprod.util import simprodtray, arguments
from icecube.simprod.util.simprodtray import RunI3Tray
from icecube import icetray, dataclasses, simclasses, sim_services
from icecube import polyplopia
from icecube import ppc
from icecube import dataio, phys_services, interfaces
from icecube.simprod import segments
from icecube.production_histograms import ProductionHistogramModule
from icecube.production_histograms.histogram_modules.simulation.mcpe_module import I3MCPEModule


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
    arguments.add_usegpus(parser, True)

    parser.add_argument("--no-RunMPHitFilter", dest="runmphitfilter",
                        default=True, action="store_false", required=False,
                        help="Don't run polyplopia's mphitfilter")
    parser.add_argument("--gpulib", dest="gpulib",
                        default="opencl", type=str, required=False,
                        help="set gpu library to load (defaults to cuda)")
    parser.add_argument("--no-volumecyl", dest="volumecyl",
                        default=True, action="store_false", required=False,
                        help="Don't set volume to regular cylinder (set flag for 300m spacing from the DOMs)")
    parser.add_argument("--MCTreeName", dest="mctreename",
                        default="I3MCTree", type=str, required=False,
                        help="Name of MCTree frame object")
    parser.add_argument("--KeepEmptyEvents", dest="keepemptyevents",
                        default=False, action="store_true", required=False,
                        help="Don't discard events with no MCPEs")
    parser.add_argument("--TempDir", dest="tempdir",
                        default=None, type=str, required=False,
                        help='Temporary working directory with the ice model')


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
    if params['propagatemuons']:
        if params['usegslrng']:
            randomServiceForPropagators = phys_services.I3GSLRandomService(seed=params['seed'] * params['nproc'] + params['procnum'])
        else:
            randomServiceForPropagators = phys_services.I3SPRNGRandomService(seed=params['seed'], nstreams=params['nproc'] * 2, streamnum=params['nproc'] + params['procnum'])

        tray.context['I3PropagatorRandomService'] = randomServiceForPropagators
        tray.AddModule("Rename", "rename_corsika_mctree", Keys=['I3MCTree', 'I3MCTree_preMuonProp'])
        tray.AddSegment(segments.PropagateMuons, 'propagator',
                        RandomService=randomServiceForPropagators,
                        **params['proposalparams'])

    tray.AddSegment(segments.PPCTraySegment, "ppc_photons",
                    GPU=params['gpu'],
                    UseGPUs=params['usegpus'],
                    UnshadowedFraction=params['efficiency'],
                    DOMOversizeFactor=params['oversize'],
                    IceModelLocation=params['icemodellocation'],
                    HoleIceParameterization=params['holeiceparametrization'],
                    IceModel=params['icemodel'],
                    volumecyl=params['volumecyl'],
                    gpulib=params['gpulib'],
                    InputMCTree=params['mctreename'],
                    keep_empty_events=params['keepemptyevents'],
                    MCPESeriesName=params['photonseriesname'],
                    tempdir=params['tempdir'])

    if params['runmphitfilter']:
        tray.AddModule("MPHitFilter", "hitfilter",
                       HitOMThreshold=1,
                       RemoveBackgroundOnly=False,
                       I3MCPESeriesMapName=params['photonseriesname'])

    if params['enablehistogram'] and params['histogramfilename']:
        tray.AddModule(ProductionHistogramModule,
                       Histograms=[I3MCPEModule],
                       OutputFilename=params['histogramfilename'])


def main():
    """
     GPU Photon propagation
    """
    # Get Params
    parser = argparse.ArgumentParser(description="PPC script")
    add_args(parser)
    params = vars(parser.parse_args())  # dict()

    # Process Params
    inputfilenamelist = [params['gcdfile']] + params['inputfilelist']

    # Execute Tray
    summary = RunI3Tray(params, configure_tray, "PPC",
                        summaryfile=params['summaryfile'],
                        inputfilenamelist=inputfilenamelist,
                        outputfile=params['outputfile'])


if __name__ == "__main__":
    main()
