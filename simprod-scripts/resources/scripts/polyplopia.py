#!/usr/bin/env python

"""
Add background coincidences to signal MC
"""

import argparse

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, simclasses, dataio
#from icecube.icetray import I3Frame
from icecube.simprod import segments
from icecube.simprod.util import arguments
from icecube.simprod.util.simprodtray import RunI3Tray
from icecube import PROPOSAL, phys_services
from icecube import polyplopia
from icecube.production_histograms import ProductionHistogramModule
from icecube.production_histograms.histogram_modules.simulation.mcpe_module import I3MCPEModule


def add_args(parser):
    """
    Args:
        parser (argparse.ArgumentParser): the command-line parser
    """
    arguments.add_gcdfile(parser)
    arguments.add_inputfile(parser)
    arguments.add_outputfile(parser)
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

    arguments.add_photonseriesname(parser)

    arguments.add_gpu(parser)
    arguments.add_usegpus(parser, True)

    parser.add_argument("--backgroundfile", dest="backgroundfile",
                        default='', type=str, required=False,
                        help='Background filename')
    parser.add_argument("--backgroundrate", dest="backgroundrate",
                        default=float('nan'), type=float, required=False,
                        help='Background rate (Hz) (don\'t use with corsika)')
    parser.add_argument("--mctype", dest="mctype",
                        default='corsika', type=str, required=False,
                        help='Type of primary simulation')
    parser.add_argument("--UsePPC", dest="useppc",
                        default=False, action="store_true", required=False,
                        help="Use PPC for photon propagation instead of CLSim.")
    parser.add_argument("--TimeWindow", dest="timewindow",
                        default=40. * I3Units.microsecond, type=float, required=False,
                        help='Coincidence time window')


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
    tray.AddSegment(segments.PolyplopiaPhotons, "coincifypes",
                    RandomService=tray.context['I3RandomService'],
                    mctype=params['mctype'],
                    mctree_name="I3MCTree",
                    bgfile=params['backgroundfile'],
                    GCDFile=params['gcdfile'],
                    timewindow=params['timewindow'],
                    GPU=params['gpu'],
                    UseGPUs=params['usegpus'],
                    UsePPC=params['useppc'],
                    IceModel=params['icemodel'],
                    IceModelLocation=params['icemodellocation'],
                    DOMOversizeFactor=params['oversize'],
                    HoleIceParameterization=params['holeiceparametrization'],
                    Efficiency=params['efficiency'],
                    PhotonSeriesName=params['photonseriesname'],
                    PROPOSALParams=params['proposalparams'],
                    rate=params['backgroundrate'] * I3Units.hertz)

    if params['histogramfilename']:
        tray.AddModule(ProductionHistogramModule,
                       Histograms=[I3MCPEModule],
                       OutputFilename=params['histogramfilename'])


def main():
    """
    Add background coincidences to signal MC
    """
    # Get Params
    parser = argparse.ArgumentParser(description="PolyplopiaModule script")
    add_args(parser)
    params = vars(parser.parse_args())  # dict()

    # Execute Tray
    summary = RunI3Tray(params, configure_tray, "PolyplopiaModule",
                        inputfilenamelist=[params['gcdfile'], params['inputfile']],
                        outputfile=params['outputfile'],
                        seed=params['seed'],
                        nstreams=params['nproc'],
                        streamnum=params['procnum'],
                        usegslrng=params['usegslrng'])


if __name__ == "__main__":
    main()
