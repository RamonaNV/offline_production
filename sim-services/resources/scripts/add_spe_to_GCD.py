#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Inject SPE Templates into a GCD file')
parser.add_argument('--fit-filename', '-f', dest='fit_filename', help='Name if input JSON file.')
parser.add_argument('--gcd', '-g', dest='gcd_filename', help='Name of GCD file.')
parser.add_argument('--as-functor', dest='as_functor', default=False, action='store_true')
args = parser.parse_args()

from I3Tray import I3Tray
from icecube import icetray, dataio

tray = I3Tray()

tray.Add("I3Reader", Filename = args.gcd_filename)

if args.as_functor:
    from icecube.phys_services.spe_fit_injector import SPEFitInjector
    spe_fit_injector = SPEFitInjector(args.fit_filename)
    tray.Add(spe_fit_injector, streams = [icetray.I3Frame.Calibration])
else:
    from icecube.phys_services.spe_fit_injector import I3SPEFitInjector
    tray.Add(I3SPEFitInjector, Filename = args.fit_filename)
    
tray.Add("I3Writer", Filename = args.gcd_filename.replace('.i3','_spe.i3'))

tray.Execute()

print("Performing a sanity check...")
import math
for fn in [args.gcd_filename, args.gcd_filename.replace('.i3','_spe.i3')]:
    with dataio.I3File(fn) as f:
        frame = f.pop_frame()
        while "I3Calibration" not in frame:
            frame = f.pop_frame()

    dc = frame['I3Calibration'].dom_cal
    nan_count = 0
    for k,v in dc.items():
        if math.isnan(v.combined_spe_charge_distribution.slc_gaus_mean):
            nan_count += 1
    print("%s NaN = %d" % (fn, nan_count))
