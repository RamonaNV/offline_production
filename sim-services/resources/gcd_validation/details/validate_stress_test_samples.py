#!/usr/bin/env python

from I3Tray import *
from icecube import icetray, dataclasses, dataio, simclasses
from os.path import expandvars

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pylab

from optparse import OptionParser
parser = OptionParser()

parser.add_option("-i","--infile", dest="INFILE",                   
                  help="Input file to read.")

parser.add_option("-m","--nhits_per_DOM", type = "int",
                  dest="nhits_per_DOM", default=20,
                  help="Number of hits per DOM")

parser.add_option("-p","--plots", action="store_true",
                  dest="GENERATE_PLOTS", default = False,
                  help="Number of hits per DOM")

(options, args) = parser.parse_args()

f = dataio.I3File(options.INFILE)

infile = dataio.I3File(options.INFILE) 
status_frame = infile.pop_frame()
while not status_frame.Has('I3DetectorStatus'):
    status_frame = infile.pop_frame()
status = status_frame.Get('I3DetectorStatus')
badDOMList = list()
badDOMListSLC = list()
if "BadDomsList" in status_frame :
    print("Found a BadDomsList in the frame.")
    print("Using this one instead.")
    badDOMList = status_frame.Get("BadDomsList")
    badDOMListSLC = status_frame.Get("BadDomsListSLC")
    print("len(badDOMList) = ",len(badDOMList))
    print("len(badDOMListSLC) = ",len(badDOMListSLC))
else:
    print(status_frame)
    try :
        from icecube.BadDomList import bad_dom_list_static
        badDOMList = bad_dom_list_static.IC86_static_bad_dom_list()
    except ImportError :
        print("ERROR : BadDomsList wasn't found in the frame")
        print("and either the BadDomList doesn't exist or")
        print("there's no static_bad_dom_list.")
        sys.exit(1)

from icecube.sim_services.sim_utils.gcd_utils import get_omgeo, get_domcal, get_domstatus
omgeo = get_omgeo( dataio.I3File(options.INFILE) )
domcal = get_domcal( dataio.I3File(options.INFILE) )
domstat = get_domstatus( dataio.I3File(options.INFILE) )
goodDOMList = [omkey for omkey,g in omgeo \
               if omkey not in badDOMList and omkey.string > 0]

counter = 0
bad_doms_with_hits = list()
while f.more():
    counter += 1
    frame = f.pop_frame()    
    
    if frame.Stop != icetray.I3Frame.DAQ : continue

    print("[  Frame %d ]" % (counter))
    print(frame)

    pulsemap = frame.Get("I3MCPulseSeriesMap")
    dlmap = frame.Get("I3DOMLaunchSeriesMap")
    calwfmap = frame.Get("CalibratedWaveforms")
    rpmap = frame.Get("WavedeformPulses")

    nhits_per_DOM = options.nhits_per_DOM

    if 'NHitsPerDOM' in frame.keys():
        print('Found `NHitsPerDOM` in frame. Override options.nhits_per_DOM')
        print()
        nhits_per_DOM = int(frame['NHitsPerDOM'].value)

    # make sure this DOM is not in the bad DOM list
    for omkey, rpseries in rpmap :
       
        charge = sum([rp.charge for rp in rpseries])
        if len(rpseries) == 0 :
            print("%s : this DOM has an empty I3RecoPulseSeries" % str(omkey))
            print("     beacon baseline ATWD0a = %f" % domcal[omkey].atwd_beacon_baseline[0,0])
            print("     beacon baseline ATWD0b = %f" % domcal[omkey].atwd_beacon_baseline[0,1])
            # how do the calibrated waveforms look?
            if options.GENERATE_PLOTS:
                atwd0 = calwfmap[omkey][0]
                fig = plt.figure()
                plt.plot(range(len(atwd0.waveform)), [v/I3Units.mV for v in atwd0.waveform])
                fig.savefig("calibrated_ATWD0_%s_%s.png" % (omkey.string, omkey.om))
                plt.clf()  

                domlaunch = dlmap[omkey][0]
                fig = plt.figure()                
                pylab.plot(range(len(domlaunch.raw_atwd[0])), [v for v in domlaunch.raw_atwd[0]])
                pylab.title("N_launches = %d LC_Bit = %s" % (len(dlmap[omkey]),domlaunch.lc_bit))
                fig.savefig("launch_ATWD0_%s_%s.png" % (omkey.string, omkey.om))
                plt.clf()  
        
        # DOMs in the badDOMListSLC should have no waveforms at all
        if omkey in badDOMListSLC :
            print("%s : this DOM is in the BAD DOM List!!!" % str(omkey))
            print("  number of recopulses = ",len(rpseries))
            print("  charge = %.2f" % charge)
            print("  number of launches = ",len(dlmap[omkey]))
            print("  lc_bit = ",dlmap[omkey][0].lc_bit)
            print("  trigger_type = ",dlmap[omkey][0].trigger_type)
            print("  trigger_mode = ",dlmap[omkey][0].trigger_mode)
            if omkey not in bad_doms_with_hits:
                bad_doms_with_hits.append(omkey)

        if(charge/float(nhits_per_DOM) < 0.2 or \
           charge/float(nhits_per_DOM) > 2.0 ) :
            print("%s : what do you think about this (%f) charge and this (%f) charge ratio? " % \
                  (str(omkey),charge,charge/float(nhits_per_DOM)))

        # The BadDOMListSLC are DOMs that are off and should not contain any hits
        # The BadDOMList are DOMs that do not participate in HLC launches
        if omkey in badDOMListSLC and omkey not in badDOMList:
            # these are SLC-only DOMs
            for dl in dlmap[omkey] :
                if dl.lc_bit :
                    print("ERROR: This %s is an SLC-only DOM with LCBit set to True." % omkey)
            
    # make sure every DOM in the good DOM list has a hit
    for omkey in goodDOMList :
        if omkey not in rpmap:
            print("%s : this DOM is good but produced no hits!!!" % str(omkey))
            print("     this is an %s DOM." % str(omgeo[omkey].omtype))

            if str(omgeo[omkey].omtype) == 'Scintillator':
                print("     No PEs were created to test the Scintillators. Skip this DOM.")
                continue

            if omkey not in pulsemap :
                print("   %s : ERROR this DOM has no PMT waveform!!!" % str(omkey))
            else:
                charge = sum([pulse.charge for pulse in pulsemap[omkey]])
                print("   %s : OK this DOM has a PMT waveform with charge %f" % (str(omkey), charge))

            if omkey not in dlmap :
                print("   %s : ERROR this DOM has no DOM launches!!!" % str(omkey))
            else:
                print("   %s : OK this DOM has %s launches." % len(dlmap[omkey]))

            if omkey not in calwfmap :
                print("   %s : ERROR this DOM has no calibrated waveforms!!!" % str(omkey))
            else:
                print("   %s : OK this DOM has %d calibrated waveforms." % len(calwfmap[omkey]))

            if omkey not in domcal :
                print("   %s : this DOM has no domcal entry!!!" % str(omkey))
            else:
                print("        impedance = %f ohms" % ( (domcal[omkey].front_end_impedance)/I3Units.ohm))
            if omkey not in domstat :
                print("   %s : this DOM has no domstat entry!!!" % str(omkey))
            else:
                print("        voltage = %f V" % ( (domstat[omkey].pmt_hv)/I3Units.V))
                print("        statusATWDa = %s" % domstat[omkey].status_atwd_a)
                print("        statusATWDb = %s" % domstat[omkey].status_atwd_b)
                print("        lcWindowPost = %s ns" % domstat[omkey].lc_window_post)

            if omkey in domcal and omkey in domstat :
                print("        gain = %f " % ( dataclasses.pmt_gain(domstat[omkey],domcal[omkey]) ))
                print("        ttime = %f ns " % ( dataclasses.transit_time(domstat[omkey],domcal[omkey])/I3Units.ns ))
            

print("number of bad DOMs with hits = ",len(bad_doms_with_hits))
print("len(badDOMList) = ",len(badDOMList))
print("len(badDOMListSLC) = ",len(badDOMListSLC))

for d in bad_doms_with_hits:
    print(d)

