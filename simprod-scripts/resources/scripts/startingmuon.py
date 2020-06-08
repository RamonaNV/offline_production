#!/usr/bin/env python

"""
IceProd Module for Simple Muon Production
"""

import os,sys
from os.path import expandvars
import logging
import argparse

from I3Tray import I3Units, I3Tray
from icecube import icetray, dataio, dataclasses
from icecube.simprod.util import simprodtray, arguments
from icecube.simprod.util.simprodtray import RunI3Tray

import random
from math import pi


def add_args(parser):
    """
    Args:
        parser (argparse.ArgumentParser): the command-line parser
    """
    arguments.add_outputfile(parser)
    arguments.add_seed(parser)

    arguments.add_nevents(parser)

    parser.add_argument("--FromEnergy", dest="fromenergy",
                        default=1.*I3Units.TeV, type=float, required=False,
                        help='Minimum energy')
    parser.add_argument("--ToEnergy", dest="toenergy",
                        default=10.*I3Units.PeV, type=float, required=False,
                        help='Maximum energy')


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
    tray.AddModule("I3InfiniteSource", Stream=icetray.I3Frame.DAQ)

    random.seed(params['seed'])

    def Generator(frame, FromEnergy = 1*I3Units.TeV, ToEnergy = 1*I3Units.TeV):
        p = dataclasses.I3Particle()
        p.energy = random.uniform(FromEnergy, ToEnergy)
        p.pos = dataclasses.I3Position(0,0,0)
        
        zenith = random.uniform(0., pi)
        azimuth = random.uniform(0., 2*pi)
        p.dir = dataclasses.I3Direction(zenith, azimuth)
        p.length = 500 * I3Units.m
        p.type = dataclasses.I3Particle.ParticleType.MuMinus
        p.location_type = dataclasses.I3Particle.LocationType.InIce
        p.time = 0. * I3Units.ns
        
        tree = dataclasses.I3MCTree()
        tree.add_primary(p)
                      
        frame["I3MCTree_preMuonProp"] = tree
         
    tray.Add(Generator,
             FromEnergy=params['fromenergy'],
             ToEnergy=params['toenergy'],
             Streams=[icetray.I3Frame.DAQ])


def main():
    """
     Injects a muon at the center of the detector.
    """
    # Get Params
    parser = argparse.ArgumentParser(description="StartingMuon script")
    add_args(parser)
    params = vars(parser.parse_args())  # dict()

    # Execute Tray
    summary = RunI3Tray(params, configure_tray, "StartingMuon",
                        outputfile=params['outputfile'],
                        executionmaxcount=params['nevents'])


if __name__ == "__main__":
    main()
