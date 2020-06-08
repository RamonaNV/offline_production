#!/usr/bin/env python


# TODO: RunId
# TODO: Argparse

import argparse
from os.path import expandvars
from copy import deepcopy

from I3Tray import I3Tray

from icecube import icetray, dataclasses, dataio
from icecube.BadDomList.BadDomListModule import BadDomListModule
from icecube import vuvuzela
from icecube.gcdserver.GCDGeneratorModule import GCDGenerator

# Notes for season 2017:
# Don't use the first runs! We had a bad calibration for DOM Coxae. Chose a run >= 129580!

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--logfile', required = False, help = "Path of the log file", type = str, default = './gcd_logfile')
parser.add_argument('-r', '--run-id', required = True, help = "Chose a run id of the season you want to generate the GCD file. Ensure that the run is a good run! In addition, check for any other probles with runs and chose one with a good calibration for all DOMs!", type = int)
args = parser.parse_args()

icetray.logging.rotating_files(args.logfile)

def remove_spe_fits(frame):
    """
    Replaces the I3Calibration frame. Removes the SPE fits from the calibration data
    since it is included in the gcdserver.
    """

    cal = deepcopy(frame['I3Calibration'])
    del frame['I3Calibration']

    for omkey, i3domcal in cal.dom_cal.items():
        i3domcal.combined_spe_charge_distribution = dataclasses.SPEChargeDistribution()
        i3domcal.mean_fadc_charge = float('nan')
        i3domcal.mean_atwd_charge = float('nan')
        i3domcal.spe_disc_calib = dataclasses.LinearFit()

        cal.dom_cal[omkey] = i3domcal

    frame['I3Calibration'] = cal

tray = I3Tray()

tray.Add(GCDGenerator, 'gcdgenerator', RunId = args.run_id)

tray.Add(remove_spe_fits, Streams = [icetray.I3Frame.Calibration])

tray.Add(BadDomListModule, 'BadDoms',
         ListName = "BadDomsList",
         DisabledKeysOnly = True,
         AddGoodSlcOnlyKeys  = True,
         Simulation = True
)

tray.Add(BadDomListModule, 'BadDomsSLC',
         ListName = "BadDomsListSLC",
         DisabledKeysOnly = True,
         AddGoodSlcOnlyKeys  = False,
         Simulation = True
)

tray.Add("Inject", "InjectNoiseParams",
	     InputNoiseFile = \
         expandvars("$I3_SRC/vuvuzela/resources/data/parameters.dat")
)

tray.Add("I3Writer","writer",
        Streams = [icetray.I3Frame.TrayInfo, \
                   icetray.I3Frame.Geometry, \
                   icetray.I3Frame.Calibration, \
                   icetray.I3Frame.DetectorStatus],
        FileName = "./gcd_snapshot.i3.gz"
)

def dump(frame):
    f = open(args.logfile, 'a')
    f.write("%s\n" % str(frame))
    
tray.AddModule(dump,"dump")

tray.Execute(4)
tray.Finish()

