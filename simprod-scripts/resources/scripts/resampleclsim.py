#!/usr/bin/env python

"""
 GPU Photon propagation w corsika resampling
"""

from icecube.simprod.util import CombineHits, DrivingTime
from icecube.simprod.util import ReadI3Summary, WriteI3Summary
from icecube.simprod.util.fileutils import download,untar,isurl
from icecube.simprod.util.simprodtray import RunI3Tray
import argparse
from icecube.simprod import segments
from icecube import icetray, dataclasses, dataio, phys_services, interfaces
from icecube import corsika_reader
from I3Tray import I3Tray, I3Units
import clsim


def add_args(parser):
    """
    Args:
        parser (argparse.ArgumentParser): the command-line parser
    """
    parser.add_argument("--OverSampling", dest="oversampling",
                        default=1, type=int, required=False,
                        help='Number of times to sample each shower.')
    parser.add_argument("--CylinderHeight", dest="cylinderheight",
                        default=1200 * I3Units.meter, type=float, required=False,
                        help='height of IceCube-centered target cylinder')
    parser.add_argument("--CylinderRadius", dest="cylinderradius",
                        default=600 * I3Units.meter, type=float, required=False,
                        help='radius of IceCube-centered target cylinder')
    clsim.add_args(parser)


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
    tray.AddModule("CORSIKAResampler", "resample",
                   OverSampling=params['oversampling'],
                   CylinderHeight=params['cylinderheight'],
                   CylinderRadius=params['cylinderradius'])

    clsim.configure_tray(tray, params, stats, logger)


def main():
    """
     GPU Photon propagation w corsika resampling
    """
    # Get Params
    parser = argparse.ArgumentParser(description="ClSimResampleCorsika script")
    add_args(parser)
    params = vars(parser.parse_args())  # dict()

    # Process Params
    inputfilenamelist = [params['gcdfile']] + params['inputfilelist']

    # Execute Tray
    summary = RunI3Tray(params, configure_tray, "ClSimResampleCorsika",
                        summaryfile=params['summaryfile'],
                        inputfilenamelist=inputfilenamelist,
                        outputfile=params['outputfile'],
                        seed=params['seed'],
                        nstreams=params['nproc'],
                        streamnum=params['procnum'],
                        usegslrng=params['usegslrng'])


if __name__ == "__main__":
    main()
