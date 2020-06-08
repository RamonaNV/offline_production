#!/usr/bin/env python
'''
A sanity checker for GCD files used in simprod.
'''
import sys
import os
import math
from os.path import expandvars
from optparse import OptionParser

from math import isnan

from icecube import icetray
from icecube import dataclasses
from icecube import simclasses

from gcd_setup import gcd_extract 

from icecube.sim_services.gcd_validation.gcd_checks import all_standard_inice_strings_exist

parser = OptionParser()
parser.add_option("-i","--inputfile", dest="GCDFILENAME",help="GCD file.")
(options, args) = parser.parse_args()

if not options.GCDFILENAME :
    print("You must specify a GCD file. (e.g. '-i <GCD_FILENAME>')")
    sys.exit(1)


gcd = gcd_extract(options.GCDFILENAME)

dom_geo_map = gcd['dom_geo_map']
dom_cal_map = gcd['dom_cal_map']
dom_status_map = gcd['dom_status_map']
vem_cal_map = gcd['vem_cal_map']
station_geo_map = gcd['station_geo_map']
bad_dom_list = gcd['bad_dom_list']
high_qe_dom_list = gcd['high_qe_dom_list']

if not all_standard_inice_strings_exist(dom_geo_map):
    print("Looks like you're missing some standard strings.")

# See ticket : http://code.icecube.wisc.edu/projects/icecube/ticket/1594
# Break this up into separate functions that test one thing.  
# Add some docs.  Read the coding standards.
# There's a great one like "Give one entity one cohesive responsibility."
for string in range(1,87):
    found_cal = False
    found_stat = False
    found_vemcal = False
    for omkey in [icetray.OMKey(string,om) for om in range(65)] :
        if omkey in dom_cal_map :
            found_cal = True
        if omkey in dom_status_map :
            found_stat = True
        if omkey.om > 60 and string < 82:
            if omkey.string not in station_geo_map:
                print('%s is missing from stationgeo' % omkey.string)
            else:
                station = station_geo_map[omkey.string]
                found_tank = False
                for tank in station:
                    if omkey in tank.omkey_list:
                        found_tank = True
                        break
                if not found_tank:
                    print('%s is missing in tanks' % omkey)
            if omkey in vem_cal_map:
                found_vemcal = True
            else :
                print('%s is missing from the VEMCal' % omkey)
        if found_cal and found_stat : continue
    if not found_cal :
        print('string %s is missing from the calibration' % string)
    if not found_stat :
        print('string %s is missing from the detector status' % string)


from icecube.vuvuzela.gcd_test import vuvuzela_test
from icecube.DOMLauncher.gcd_test import pmt_response_sim_test
from icecube.DOMLauncher.gcd_test import dom_launcher_test
from icecube.topsimulator.gcd_test import topsimulator_test
from icecube.trigger_sim.gcd_test import trigger_sim_test

def photon_propagators_test(omkey, i3omgeo, domcal, domstat):
    if i3omgeo.omtype != dataclasses.I3OMGeo.OMType.IceCube :
        return True
    
    if isnan(i3omgeo.position.x) \
       or isnan(i3omgeo.position.y) \
       or isnan(i3omgeo.position.z) \
       or isnan(domcal.relative_dom_eff) \
       or domstat.pmt_hv <= 0. :

        print('  %s  photon propagation !!' % (str(omkey)) )
        print('     pos = %.2f' % i3omgeo.position )
        print('     RDE = %.2f' % domcal.relative_dom_eff )
        print('     PMT HV = %.2f' % domstat.pmt_hv )

        return False
    return True


# This is True if every DOM passes individual tests.
all_pass = True
for omkey, i3omgeo in dom_geo_map:

    if omkey not in bad_dom_list \
           and omkey in dom_cal_map \
           and omkey in dom_status_map:
        print("Testing DOM %s" % omkey)
        domcal = dom_cal_map[omkey]
        domstat = dom_status_map[omkey]

        pass_vuvuzela = vuvuzela_test(omkey, i3omgeo, domcal)
        pass_pmt = pmt_response_sim_test(omkey, domcal, domstat)
        pass_dom_launcher = dom_launcher_test(omkey, i3omgeo, domcal, domstat)
        pass_photon_propagators = photon_propagators_test(omkey, i3omgeo, domcal, domstat)

        if i3omgeo.omtype == dataclasses.I3OMGeo.OMType.IceTop:
            if omkey.string in station_geo_map:
                station = station_geo_map[omkey.string]
                for tank in station:
                    if omkey in tank.omkey_list:
                        vemcal = vem_cal_map[omkey]
                        pass_top_sim = topsimulator_test(omkey, tank, domcal, vemcal, domstat)
                        break
        else:
            pass_top_sim = True

        if not pass_vuvuzela : 
            print ("FAIL : Vuvuzela")
        if not pass_pmt : 
            print ("FAIL : PMTResponseSimulator")
        if not pass_dom_launcher : 
            print ("FAIL : DOMLauncher")
        if not pass_photon_propagators : 
            print ("FAIL : ppc and clsim")
        if not pass_top_sim : 
            print ("FAIL : I3TopSimulator")

        all_pass = all_pass \
            and pass_vuvuzela \
            and pass_pmt \
            and pass_photon_propagators \
            and pass_dom_launcher \
            and pass_top_sim

# report back to the mothership
if not all_pass :
    sys.exit(1)
else :
    sys.exit(0)        

