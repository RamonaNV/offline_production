#!/usr/bin/env python

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataio, dataclasses, phys_services
from icecube.simprod.util import ReadI3Summary, WriteI3Summary, BasicCounter, DAQCounter
from icecube.simprod.util import simprodtray, arguments
from icecube.simprod.util.simprodtray import RunI3Tray
import argparse
from icecube import neutrino_generator, earthmodel_service, PROPOSAL, cmc
from icecube.simprod import segments
import json
from icecube.production_histograms import ProductionHistogramModule
from icecube.production_histograms.histogram_modules.simulation.mctree_primary import I3MCTreePrimaryModule
from icecube.production_histograms.histogram_modules.simulation.mctree import I3MCTreeModule


def add_args(parser):
    """
    Args:
        parser (argparse.ArgumentParser): the command-line parser
    """
    arguments.add_gcdfile(parser,required=False)
    arguments.add_outputfile(parser)
    arguments.add_summaryfile(parser)
    arguments.add_enablehistogram(parser)
    arguments.add_histogramfilename(parser)

    arguments.add_nproc(parser)
    arguments.add_procnum(parser)
    arguments.add_seed(parser)
    arguments.add_usegslrng(parser)

    arguments.add_nevents(parser)

    arguments.add_proposalparams(parser)
    arguments.add_propagatemuons(parser, False)

    parser.add_argument("--runid", dest="runid",
                        default=None, type=int, required=False,
                        help='Run number to use in S-Frame')    
    parser.add_argument("--SimMode", dest="simmode",
                        default='FULL', type=str, required=False,
                        help='simulation mode')
    parser.add_argument("--VTXGenMode", dest="vtxgenmode",
                        default='NuGen', type=str, required=False,
                        help='currently supports only NuGen')
    parser.add_argument("--InjectionMode", dest="injectionmode",
                        default='Surface', type=str, required=False,
                        help='injection mode')
    parser.add_argument("--CylinderParams", dest="cylinderparams",
                        default=[0, 0, 0, 0, 0], type=arguments.float_comma_list, required=False,
                        help='For CIRCLE[radius, active_height_before, active_height_after],'' for SURFACE[radius, length, center_x, center_y, center_z]')
    parser.add_argument("--no-AutoExtendMuonVolume", dest="autoextendmuonvolume",
                        default=True, action="store_false", required=False,
                        help="Don't use detection volume extension (set flag for starting events)")
    parser.add_argument("--NuFlavor", dest="nuflavor",
                        default='', type=str, required=False,
                        help='Use Legacy injector : Neutrino Flavor')
    parser.add_argument("--NuTypes", dest="nutypes",
                        default=dict(), type=json.loads, required=False,
                        help='Use new injector. Set dictionary of neutrino types and ratio. e.g. {NuMu:1, NuMuBar:1, NuEBar:1}')
    parser.add_argument("--Polyplopia", dest="polyplopia",
                        default=False, action="store_true", required=False,
                        help='Produce coincident showers')
    parser.add_argument("--gamma", dest="gamma",
                        default=2.0, type=float, required=False,
                        help='Gamma index')
    parser.add_argument("--FromEnergy", dest="fromenergy",
                        default=1. * I3Units.TeV, type=float, required=False,
                        help='Minimum energy')
    parser.add_argument("--ToEnergy", dest="toenergy",
                        default=10. * I3Units.PeV, type=float, required=False,
                        help='Maximum energy')
    parser.add_argument("--ZenithMin", dest="zenithmin",
                        default=0., type=float, required=False,
                        help='min zenith')
    parser.add_argument("--ZenithMax", dest="zenithmax",
                        default=180. * I3Units.degree, type=float, required=False,
                        help='max zenith')
    parser.add_argument("--AzimuthMin", dest="azimuthmin",
                        default=0., type=float, required=False,
                        help='min azimuth')
    parser.add_argument("--AzimuthMax", dest="azimuthmax",
                        default=360. * I3Units.degree, type=float, required=False,
                        help='max azimuth')
    parser.add_argument("--ZenithSamplingMode", dest="zenithsamplingmode",
                        default='ANGEMU', type=str, required=False,
                        help='zenith sampling mode')
    parser.add_argument("--UseDifferentialXsection", dest="usedifferentialxsection",
                        default=False, action="store_true", required=False,
                        help='Do you use differential cross sections?')
    parser.add_argument("--CrossSections", dest="crosssections",
                        default='csms', type=str, required=False,
                        help='cross section files')
    parser.add_argument("--CrossSectionsPath", dest="crosssectionspath",
                        default='', type=str, required=False,
                        help='cross section tables path')
    parser.add_argument("--ParamsMap", dest="paramsmap",
                        default=dict(), type=json.loads, required=False,
                        help='any other parameters')
    parser.add_argument("--BackgroundFile", dest="backgroundfile",
                        default="", type=str, required=False,
                        help='pre-generated coincident showers file')


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
    if params['usegslrng']:
        randomServiceForPropagators = phys_services.I3GSLRandomService(seed=2 * (params['seed'] * params['nproc'] + params['procnum']))
    else:
        randomServiceForPropagators = phys_services.I3SPRNGRandomService(seed=params['seed'], nstreams=params['nproc'] * 2, streamnum=params['nproc'] + params['procnum'])

    tray.AddModule("I3InfiniteSource", "TheSource",
                   Prefix=params['gcdfile'],
                   Stream=icetray.I3Frame.DAQ)
    tray.AddModule(DAQCounter, "counter3", nevents=int(params['nevents']))

    vnutypes = []
    vtyperatio = []
    for key in params['nutypes']:
        vnutypes.append(key)
        vtyperatio.append(params['nutypes[key]'])

    tray.AddSegment(segments.GenerateNeutrinos, 'generator',
                    RandomService=tray.context['I3RandomService'],
                    RunID=params['runid'],
                    NumEvents=params['nevents'],
                    SimMode=params['simmode'],
                    VTXGenMode=params['vtxgenmode'],
                    InjectionMode=params['injectionmode'],
                    CylinderParams=params['cylinderparams'],
                    AutoExtendMuonVolume=params['autoextendmuonvolume'],
                    Flavor=params['nuflavor'],
                    NuTypes=vnutypes,
                    PrimaryTypeRatio=vtyperatio,
                    GammaIndex=params['gamma'],
                    FromEnergy=params['fromenergy'],
                    ToEnergy=params['toenergy'],
                    ZenithRange=[params['zenithmin'], params['zenithmax']],
                    AzimuthRange=[params['azimuthmin'], params['azimuthmax']],
                    UseDifferentialXsection=params['usedifferentialxsection'],
                    CrossSections=params['crosssections'],
                    CrossSectionsPath=params['crosssectionspath'],
                    ZenithSamplingMode=params['zenithsamplingmode'],
                    ParamsMap=params['paramsmap'])

    if params['polyplopia']:
        tray.AddSegment(segments.PolyplopiaSegment, "coincify",
                        RandomService=tray.context['I3RandomService'],
                        mctype='NuGen',
                        bgfile=params['backgroundfile'],
                        timewindow=40. * I3Units.microsecond,
                        rate=5.0 * I3Units.kilohertz)

    if params['propagatemuons']:
        tray.context['I3PropagatorRandomService'] = randomServiceForPropagators
        tray.AddSegment(segments.PropagateMuons, 'propagator',
                        RandomService=randomServiceForPropagators,
                        **params['proposalparams'])

    tray.AddModule(BasicCounter, "count_g",
                   Streams=[icetray.I3Frame.DAQ],
                   name="Generated Events",
                   Stats=stats)

    if params['enablehistogram'] and params['histogramfilename']:
        tray.AddModule(ProductionHistogramModule,
                       Histograms=[I3MCTreePrimaryModule, I3MCTreeModule],
                       OutputFilename=params['histogramfilename'])


def main():
    # Get Params
    parser = argparse.ArgumentParser(description="NuGen script")
    add_args(parser)
    params = vars(parser.parse_args())  # dict()

    # Execute Tray
    summary = RunI3Tray(params, configure_tray, "NuGen",
                        summaryfile=params['summaryfile'],
                        summaryin=dataclasses.I3MapStringDouble(),
                        outputfile=params['outputfile'],
                        seed=params['seed'],
                        nstreams=params['nproc'],
                        streamnum=params['procnum'],
                        usegslrng=params['usegslrng'])


if __name__ == "__main__":
    main()
