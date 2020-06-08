#!/usr/bin/env python

"""
Common argparse arguments and utility functions
"""

from os.path import expandvars
import json


# split up a comma-delimited string into list
def _comma_list(arg_list, type):
    return [type(a) for a in arg_list.split(',')]

def int_comma_list(arg_list):
    return _comma_list(arg_list, int)

def float_comma_list(arg_list):
    return _comma_list(arg_list, float)

def str_comma_list(arg_list):
    return _comma_list(arg_list, str)
    

# common arguments    
def add_gcdfile(parser,required=False):
    parser.add_argument("--gcdfile", dest="gcdfile",
                        default='', type=str, required=required,
                        help='GeoCalibDetStatus filename')
                        
def add_inputfile(parser):
    parser.add_argument("--inputfile", dest="inputfile",
                        default='', type=str, required=True,
                        help='Input filename')
               
def add_outputfile(parser):         
    parser.add_argument("--outputfile", dest="outputfile",
                        default='', type=str, required=True,
                        help='Output filename')

def add_inputfilelist(parser):
    parser.add_argument("--inputfilelist", dest="inputfilelist",
                        default=[], type=str_comma_list, required=True,
                        help='List of input filenames')

def add_nproc(parser):
    parser.add_argument("--nproc", dest="nproc",
                        default=1, type=int, required=False,
                        help='Number of jobs (Number of RNG streams)')

def add_procnum(parser):
    parser.add_argument("--procnum", dest="procnum",
                        default=0, type=int, required=False,
                        help='job number (RNG stream number)')

def add_seed(parser):
    parser.add_argument("--seed", dest="seed",
                        default=0, type=int, required=False,
                        help='RNG seed')

def add_usegslrng(parser):
    parser.add_argument("--UseGSLRNG", dest="usegslrng",
                        default=False, action="store_true", required=False,
                        help="Use I3GSLRandomService")

def add_enablehistogram(parser):
    parser.add_argument("--EnableHistogram", dest="enablehistogram",
                        default=False, action="store_true", required=False,
                        help='Write a SanityChecker histogram file.')

def add_histogramfilename(parser):
    parser.add_argument("--HistogramFilename", dest="histogramfilename",
                        default=None, type=str, required=False,
                        help='Histogram filename.')

def add_summaryfile(parser):
    parser.add_argument("--summaryfile", dest="summaryfile",
                        default='summary.json', type=str, required=False,
                        help='JSON Summary filename')

def add_nevents(parser):
    parser.add_argument("--nevents", dest="nevents",
                        default=1, type=int, required=False,
                        help='Number of events')

def add_gpu(parser):
    parser.add_argument("--GPU", dest="gpu",
                        default=None, type=str, required=False,
                        help="Graphics Processing Unit number (shoud default to environment if None)") 

def add_icemodellocation(parser):
    parser.add_argument("--IceModelLocation", dest="icemodellocation",
                        default=expandvars("$I3_BUILD/ice-models/resources/models"),
                        type=str, required=False,
                        help="Location of ice model param files")

def add_icemodel(parser):
    parser.add_argument("--IceModel", dest="icemodel",
                        default="spice_3.2", type=str, required=False,
                        help="ice model subdirectory")

def add_proposalparams(parser):
    parser.add_argument("--PROPOSALParams", dest="proposalparams",
                        default=dict(), type=json.loads, required=False,
                        help='any other parameters for proposal')

def add_holeiceparametrization(parser):
    parser.add_argument("--holeiceparametrization", dest="holeiceparametrization",
                        default=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm"),
                        type=str, required=False,
                        help="Location of hole ice param files")

def add_oversize(parser):
    parser.add_argument("--oversize", dest="oversize",
                        default=5, type=int, required=False,
                        help="over-R: DOM radius oversize scaling factor")

def add_efficiency(parser):
    parser.add_argument("--efficiency", dest="efficiency",
                        default=1.00, type=float, required=False,
                        help="overall DOM efficiency correction")

def add_photonseriesname(parser):
    parser.add_argument("--PhotonSeriesName", dest="photonseriesname",
                        default="I3MCPESeriesMap", type=str, required=False,
                        help="Photon Series Name")

def add_propagatemuons(parser, default):
    if default:
        parser.add_argument("--no-PropagateMuons", dest="propagatemuons",
                            default=True, action="store_false", required=False,
                            help="Don't run PROPOSAL to do in-ice propagation")
    else:
        parser.add_argument("--PropagateMuons", dest="propagatemuons",
                            default=False, action="store_true", required=False,
                            help="Run PROPOSAL to do in-ice propagation")

def add_usegpus(parser, default):
    if default:
        parser.add_argument("--no-UseGPUs", dest="usegpus",
                            default=True, action="store_false", required=False,
                            help="Don't use Graphics Processing Unit for photon propagation")
    else:
        parser.add_argument("--UseGPUs", dest="usegpus",
                            default=False, action="store_true", required=False,
                            help="Use Graphics Processing Unit for photon propagation")


