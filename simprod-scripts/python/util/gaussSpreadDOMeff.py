#!/usr/bin/env python
import sys, os
from os.path import expandvars
from math import isnan
import random
from .. import ipmodule

class JitterGCD(ipmodule.IPBaseClass):
   def __init__(self):

      ipmodule.IPBaseClass.__init__(self)
      self.AddParameter('inputfile','input GCD','')
      self.AddParameter('outputfile','input GCD','')


   def Execute(self,stats):
      from I3Tray import I3Tray, I3Units
      from icecube import icetray, dataclasses, dataio, simclasses
      
      from icecube.BadDomList import bad_dom_list_static
      badOMs = bad_dom_list_static.IC79_static_bad_dom_list()

      infile  = self.GetParameter('inputfile')
      outfile = self.GetParameter('outputfile')

      infile  = dataio.I3File(infile)
      geo_frame = infile.pop_frame()
      while not geo_frame.Has('I3Geometry'): geo_frame = infile.pop_frame()
      geometry = geo_frame.Get('I3Geometry')

      cal_frame = infile.pop_frame()
      while not cal_frame.Has('I3Calibration'): cal_frame = infile.pop_frame()
      calibration = cal_frame.Get('I3Calibration')

      status_frame = infile.pop_frame()
      while not status_frame.Has('I3DetectorStatus'): status_frame = infile.pop_frame()
      status = status_frame.Get('I3DetectorStatus')

      dom_geo = geometry.omgeo
      dom_cal = calibration.domCal  #.dom_cal
      dom_status = status.domStatus #.dom_status
	
      c_and_d_strings_to_check = range(1,87)
      low_noise_DOMs_l = [ icetray.OMKey(82,54),  icetray.OMKey(84,54),  icetray.OMKey(85,55)]

      strings_IC86 = [ 1, 7, 14, 22, 31, 79, 80 ]
      strings_IC79 = range(1,87)

      for e,p in dom_geo:
          if e not in badOMs and e in dom_cal and e in dom_status:				
               cal_this_om = dom_cal[e]
               status_this_om = dom_status[e]

          #### Numbers are taken from Henrik Johanson's licentiate thesis (2008). 
          #    "Calibration of photon detection efficiency in IceCube"
          #    from data he concludes for the DOM RDE: Gaussian distribution with central value 1, and sigma=0.08673

          if (cal_this_om.RelativeDomEff>1.3) :	#.relative_dom_eff
             check = cal_this_om.RelativeDomEff
             rand = random.gauss(1.35, 0.08673)
             calibration.domCal[e].RelativeDomEff = rand
             icetray.logging.log_debug("  correcting RDE from %.2f to %.2f in %s" %
                                       (check,calibration.domCal[e].RelativeDomEff,e))

          if (cal_this_om.RelativeDomEff==1.0) :
             check = cal_this_om.RelativeDomEff
             rand = random.gauss(1.0, 0.08673)
             calibration.domCal[e].RelativeDomEff = rand
             icetray.logging.log_debug("  correcting RDE from %.2f to %.2f in %s" %
                                       (check,calibration.domCal[e].RelativeDomEff,e))
             
      del geo_frame['I3Geometry']
      geo_frame['I3Geometry'] = geometry

      del cal_frame['I3Calibration']
      cal_frame['I3Calibration'] = calibration

      del status_frame['I3DetectorStatus']
      status_frame['I3DetectorStatus'] = status


      infile.close()
      outfile = dataio.I3File(outfile, dataio.I3File.Mode.Writing)
      outfile.push(geo_frame)
      outfile.push(cal_frame)
      outfile.push(status_frame)
      outfile.close()

      return 0
