import time
script_start_time = time.time()

import time, os, sys, json
import numpy as np
import argparse

#:--
# PREPARATIONS
#:--

# Argument parsing below

parser = argparse.ArgumentParser(description = "Simulating monopoles for the MMACT analysis")

parser.add_argument( "-n", "--numevents", dest="nev",       type=int, default=-1, help="The number of events to simulate."                )
parser.add_argument( "-d", "--outdir",    dest="outdir",    type=str, default="", help="The directory where the files should be created." )
parser.add_argument( "-p", "--iprocess",  dest="iprocess",  type=int, default=-1, help="The number of the current process."               )
parser.add_argument( "-r", "--nprocess",  dest="nprocess",  type=int, default=-1, help="The total number of processes."                   )
parser.add_argument( "-g", "--gcdfile",   dest="gcdfile",   type=str, default="", help="The path to the GCD file."                        )
parser.add_argument( "-m", "--icemodel",  dest="icemodel",  type=str, default="", help="The path to the icemodel."                        )
parser.add_argument( "-v", "--verbose",   dest="verbose",   type=int, default=1,  help="How verbose is this script?"                      )

args = parser.parse_args()

import mmact_utils as utils
mmp = utils.mmact_print(script_start_time,args.verbose)
mmp.start()
mmp.vbprint( "Parsed Arguments",  [1,2] )


from icecube import icetray, dataclasses, dataio, sim_services, phys_services
from icecube import monopole_generator
from icecube.icetray import I3Units
from icecube.simprod.segments import DetectorSim, Calibration, PropagatePhotons, PropagateMuons, PPCTraySegment

import I3Tray
from I3Tray import *
from I3Tray import load


mmp.label = "MMACT"
mctreename = "I3MCTree"

mmp.vbprint( "Imported Stuff",    [1,2] )

#load("xppc")
#load("ppc")
#load("libDOMLauncher")

mmp.vbprint( "Loaded Stuff",      [1,2] )

#exit("Change the meta-project in the shebang line.")



#:--
# SETTING UP SOME PARAMETERS
#:--

temp_filename = utils.filename_template
temp_filename = temp_filename.replace( "FLAVOR",        "monopole"                                                      )
temp_filename = temp_filename.replace( "BETALOW",       "{:04d}".format(int(utils.gen_params["beta_spectrum"][0]*1000)) )
temp_filename = temp_filename.replace( "BETAHIGH",      "{:04d}".format(int(utils.gen_params["beta_spectrum"][1]*1000)) )
temp_filename = temp_filename.replace( "PROCESSNUMBER", "{:04d}".format(int(args.iprocess))                                  )
outname_gen   = temp_filename.replace( "DATALEVEL",     "generator"                                                     )
outname_trigg = temp_filename.replace( "DATALEVEL",     "trigger"                                                       )

gcdname  = utils.default_settings["gcd"]      if not args.gcdfile  else args.gcdfile
icemodel = utils.default_settings["icemodel"] if not args.icemodel else args.icemodel
nev      = utils.default_settings["n_events"] if  0 >= args.nev      else args.nev

if args.iprocess<0 or args.nprocess<0 or args.iprocess>=args.nprocess:
	exit("At least one invalid processing number input!")

#os.putenv("PPCTABLESDIR",os.path.expandvars(icemodel))

#:--
# STARTING THE TRAY
#:--

tray = I3Tray()

mmp.vbprint( "Started the tray.", [1,2] )

randomService = phys_services.I3SPRNGRandomService(
             		seed = 0,
             		nstreams = args.nprocess,
             		streamnum = args.iprocess)

tray.context['I3RandomService'] = randomService

tray.AddModule( "I3InfiniteSource", Prefix = gcdname )

mmp.vbprint( "Added random service to the tray.", [1,2] )

# Generation below

tray.AddModule("I3MonopoleGenerator",
	Mass       = utils.gen_params["monopole_mass"],
	BetaRange  = utils.gen_params["beta_spectrum"],
	Disk_dist  = utils.gen_params["disk_distance"],
	Disk_rad   = utils.gen_params["disk_radius"],
	)
mmp.vbprint( "Added monopole generator to the tray.", [1,2] )

# Propagation below

tray.AddModule("I3MonopolePropagator",
	MaxDistanceFromCenter = utils.gen_params["dist_to_cent_max"],
	MaxLength             = utils.gen_params["step_length_max"],
	StepSize              = np.nan,
	)
mmp.vbprint( "Added monopole propagator to the tray.", [1,2] )

tray.Add(utils.check_monopole_lengths_10m, streams = [icetray.I3Frame.DAQ])
mmp.vbprint( "Added check_monopole_lengths_10m to the tray.", [1,2] )

tray.Add(utils.add_MCPrimaryParticle, streams = [icetray.I3Frame.DAQ])
mmp.vbprint( "Added add_MCPrimaryParticle to the tray.", [1,2] )

randomServiceForPropagators = phys_services.I3SPRNGRandomService(
             		seed = 0,
             		nstreams = args.nprocess*2,
             		streamnum = args.nprocess + args.iprocess)
tray.context['I3PropagatorRandomService'] = randomServiceForPropagators

tray.AddModule("Rename","rename_corsika_mctree",
        	           Keys=[mctreename,mctreename+'_preMuonProp'])
tray.AddSegment(PropagateMuons, 'propagator',
                        RandomService= randomServiceForPropagators
) 


# Light production below
tray.AddSegment(PPCTraySegment, "makeCLSimHits",
            UseGPUs = True,
            IceModelLocation = os.path.expandvars("$I3_SRC/ice-models/resources/models"),
            IceModel = icemodel,
            UnshadowedFraction = 1.,
            DOMOversizeFactor = 1.,
            HoleIceParameterization = os.path.expandvars("$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm"),
            InputMCTree="I3MCTree",
            MCPESeriesName = "I3MCPESeriesMap")
mmp.vbprint( "Added PPC to the tray.", [1,2] )

# Unify name for rest of icecube

tray.Add( "Rename", "PPCRename", Keys=["MCPESeriesMap", "I3MCPESeriesMap"] )

tray.Add( DetectorSim, "DetectorTrigger",
		  RandomService        = 'I3RandomService',
		  GCDFile              = args.gcdfile,
		  InputPESeriesMapName = 'I3MCPESeriesMap',
		  SkipNoiseGenerator   = False,
		  RunID                = 0,
		  KeepPropagatedMCTree = True,
		  KeepMCHits           = False,
		  KeepMCPulses         = False,
		  LowMem               = False,
		  BeaconLaunches       = True,
		  TimeShiftSkipKeys    = [],
		  FilterTrigger        = True
	)
mmp.vbprint( "Added DetectorSim to the tray.", [1,2] )


# Writing out the files

tray.AddModule( "I3Writer",
	filename = outname_trigg,
	streams  = [ icetray.I3Frame.DAQ, icetray.I3Frame.Physics] ,
	SkipKeys = [ "I3MCPulseSeriesMapPrimaryIDMap" ],
	)
mmp.vbprint( "Added I3Writer to the tray, to write out the trigger level I3 file.", [1,2] )


# Finishing and executing

tray.AddModule("TrashCan","Trash") # perhaps not needed anymore
mmp.vbprint( "TrashCan added.", [1,2] )
mmp.vbprint( "Execution upcoming!", [1,2] )
tray.Execute( nev+3 ) # Execute(n_events_you_want+4), the four are I, G, C, D (trayinfo, geometry, calibration, detector status)
mmp.vbprint( "Execution done!", [1,2] )
tray.Finish()
mmp.vbprint( "Tray finished!", [1,2] )
del tray
mmp.vbprint( "Tray deleted!", [1,2] )

mmp.finish()
