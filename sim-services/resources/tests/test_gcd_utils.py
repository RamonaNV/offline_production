#!/usr/bin/env python

import os
import unittest
from icecube import dataio
from icecube import dataclasses
from icecube.sim_services.sim_utils import gcd_utils


class TestGCDUtils(unittest.TestCase):

    def setUp(self):
        filename = os.path.expandvars("$I3_TESTDATA/GCD/GeoCalibDetectorStatus_2013.56429_V1.i3.gz")
        self.gcd_file = dataio.I3File(filename)
    
    def test_get_time(self):
        obj = gcd_utils.get_omgeo(self.gcd_file)
        self.assertNotEqual(obj, None, "fail")
                
    def test_get_omgeo(self):
        obj = gcd_utils.get_domcal(self.gcd_file)
        self.assertNotEqual(obj, None, "fail")

    def test_get_domcal(self):
        obj = gcd_utils.get_domcal(self.gcd_file)
        self.assertNotEqual(obj, None, "fail")

    def test_get_domstatus(self):
        obj = gcd_utils.get_domstatus(self.gcd_file)
        self.assertNotEqual(obj, None, "fail")

    def test_get_triggerstatus(self):
        obj = gcd_utils.get_triggerstatus(self.gcd_file)
        self.assertNotEqual(obj, None, "fail")

    def test_put_triggerstatus(self):
        tkey = dataclasses.TriggerKey()
        tkey.source = dataclasses.IN_ICE
        tkey.type = dataclasses.SIMPLE_MULTIPLICITY
        tkey.config_id = 999

        tstat = dataclasses.I3TriggerStatus()
        tstat.trigger_name = "MyNewTrigger"
        tstat.trigger_settings["threshold"] = '3'
        tstat.trigger_settings["timeWindow"] = '2500'
        tstat.trigger_settings["domSet"] = '5'

        roc = dataclasses.I3TriggerStatus.I3TriggerReadoutConfig()
        roc.readout_time_minus = 10000
        roc.readout_time_plus = 10000
        roc.readout_time_offset = 0

        # add the readout config to the trigger status
        tstat.readout_settings[dataclasses.I3TriggerStatus.ALL] = roc

        tstat_map = dataclasses.I3TriggerStatusMap()
        tstat_map[tkey] = tstat
        
        obj = gcd_utils.put_triggerstatus(tstat_map, self.gcd_file, "./GCD_test_utils.i3.gz")
        self.assertNotEqual(obj, None, "fail")    
        
unittest.main()
