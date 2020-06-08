#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--infile",
                  dest="INFILE", default = None,
                  help="GCD file to correct.")

parser.add_option("-o","--outfile",
                  dest="OUTFILE", default = None,
                  help="Corrected GCD file.")

parser.add_option("-l","--logfile",
                  dest="LOGFILE", default="./gcd_logfile" ,
                  help="Name of logfile.")

(options, args) = parser.parse_args()

import os
import sys
from os.path import expandvars
from math import isnan

from I3Tray import I3Units
from icecube import icetray
from icecube import dataclasses
from icecube import dataio
from icecube import simclasses
from icecube.icetray import OMKey

logfile = open(options.LOGFILE, "a")

if not options.INFILE :
    msg = "ERROR : You must specify and input file."
    logfile.write(msg)
    print(msg)
    sys.exit(1)

if not options.OUTFILE :
    msg = "ERROR : You must specify and output file."
    logfile.write(msg)
    print(msg)    
    sys.exit(1)

infile = dataio.I3File(options.INFILE)

geo_frame = infile.pop_frame()
while not geo_frame.Has('I3Geometry') :
    geo_frame = infile.pop_frame()
geometry = geo_frame.Get('I3Geometry')

cal_frame = infile.pop_frame()
while not cal_frame.Has('I3Calibration') :
    cal_frame = infile.pop_frame()
calibration = cal_frame.Get('I3Calibration')

status_frame = infile.pop_frame()
while not status_frame.Has('I3DetectorStatus') :
    status_frame = infile.pop_frame()
status = status_frame.Get('I3DetectorStatus')

badOMs = list()
if "BadDomsList" in status_frame :
    logfile.write("Found a BadDomsList in the frame.  Gonna use it.\n")
    badOMs = status_frame.Get("BadDomsList")
else:
    print(status_frame)
    try :
        from icecube.BadDomList import bad_dom_list_static
        badOMs = bad_dom_list_static.IC86_static_bad_dom_list()
    except ImportError :        
        logfile.write("ERROR : BadDomsList wasn't found in the frame\n")
        logfile.write("and either the BadDomList doesn't exist or\n")
        logfile.write("there's no static_bad_dom_list.\n")
        sys.exit(1)

dom_geo = geometry.omgeo
dom_cal = calibration.dom_cal
dom_status = status.dom_status

low_noise_DOMs_l = [ icetray.OMKey(82,54),  \
                     icetray.OMKey(84,54),  \
                     icetray.OMKey(85,55)]

strings_IC86 = [ 1, 7, 14, 22, 31, 79, 80 ]

high_QE = [ icetray.OMKey( 79, i ) for i in range(30, 45) \
            if i not in [ 32, 41, 43 ] ]
high_QE.extend( [ icetray.OMKey( 80, i ) for i in range(30, 44) ] )
#Mark Krasberg (post deployment after surface testing)
high_QE.extend( [ icetray.OMKey( 43, 55 ) ] )
for s in range(81,87):
    high_QE.extend( [ icetray.OMKey( s, i ) for i in range(1, 61) ] )

for omkey, omgeo in dom_geo:

    if omkey not in badOMs \
           and omkey in dom_cal \
           and omkey in dom_status:

        domcal = dom_cal[omkey]
        domstat = dom_status[omkey]

        threshold = dataclasses.spe_pmt_threshold(domstat, domcal)
        threshold /= I3Units.mV

        if threshold < 0 :
            logfile.write('Pathological PMT discriminator threshold')
            logfile.write('  %s  threshold = %2.2f mV' % (str(omkey), threshold))
            fit = dataclasses.LinearFit()
            fit.slope = NaN
            fit.intercept = NaN
            calibration.dom_cal[omkey].PMTDiscCalib = fit

            dc = calibration.dom_cal[omkey]
            thresh = dataclasses.spe_pmt_threshold(domstat, dc)
            logfile.write('  correcting to %2.2f mV' % thresh/I3Units.mV )

        if omgeo.omtype == dataclasses.I3OMGeo.IceCube :
            if omkey not in badOMs and isnan( domcal.dom_noise_rate ) :
                logfile.write("ERROR : no valid noise rate for this DOM %s " % str(omkey))
                logfile.write("  NOT CORRECTING")

            if isnan(domcal.relative_dom_eff) :
                if omkey.string in strings_IC86 :
                    new_RDE = 1.35 if e in high_QE else 1.0
                    calibration.dom_cal[omkey].relative_dom_eff = new_RDE                    
                    logfile.write(" correcting RDE from 'nan' to %.2f in %s" % \
                          (new_RDE,omkey))

            # check for unusually low noise DOMs that were
            # incorrectly translated into the DB
            if omkey in low_noise_DOMs_l :
                noiseRate = calibration.dom_cal[omkey].dom_noise_rate
                if noiseRate < 400 * I3Units.hertz :
                    new_rate = noiseRate + 1*I3Units.kilohertz
                    calibration.dom_cal[omkey].dom_noise_rate = new_rate
                    logfile.write("  correcting noise from %fHz to %fHz in %s" % \
                          (noiseRate/I3Units.hertz,
                           new_rate/I3Units.hertz,omkey))
                           


# let's clean the trigger cruft out
# we only simulate InIce and IceTop trigggers and
# only SM, Cluster, Cylinder, and SlowMonopole
# so that's all we're keeping
for tkey, ts in status.trigger_status :
    if ( tkey.source != dataclasses.IN_ICE \
         and tkey.source != dataclasses.ICE_TOP) or \
         ( tkey.type != dataclasses.SIMPLE_MULTIPLICITY and
           tkey.type != dataclasses.STRING and
           tkey.type != dataclasses.VOLUME and
           tkey.type != dataclasses.SLOW_PARTICLE ) :
        del status.trigger_status[tkey]

# now we clean out the AMANDA geometry as well
for omkey, i3omgeo in geometry.omgeo:
  if i3omgeo.omtype == dataclasses.I3OMGeo.AMANDA :
    logfile.write("Removing AMANDA OM %s from the geometry." % str(omkey))
    del geometry.omgeo[omkey]

del geo_frame['I3Geometry']
geo_frame['I3Geometry'] = geometry

del cal_frame['I3Calibration']
cal_frame['I3Calibration'] = calibration

del status_frame['I3DetectorStatus']
status_frame['I3DetectorStatus'] = status

outfile = dataio.I3File(options.OUTFILE, dataio.I3File.Mode.Writing)
outfile.push(geo_frame)
outfile.push(cal_frame)
outfile.push(status_frame)
outfile.close()

