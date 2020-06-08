#!/usr/bin/env python

from icecube import icetray
from icecube import dataio

def gcd_extract(gcd_filename):
    
    gcdfile = dataio.I3File(gcd_filename)

    frame = gcdfile.pop_frame()
    while not frame.Has('I3Geometry'):
        frame = gcdfile.pop_frame()
    geometry = frame.Get('I3Geometry')

    while not frame.Has('I3Calibration'):
        frame = gcdfile.pop_frame()
    calibration = frame.Get('I3Calibration')

    status_frame = gcdfile.pop_frame()
    while not status_frame.Has('I3DetectorStatus'):
        status_frame = gcdfile.pop_frame()
    status = status_frame.Get('I3DetectorStatus')

    dom_geo_map = geometry.omgeo
    dom_cal_map = calibration.dom_cal
    dom_status_map = status.dom_status
    trigger_status_map = status.trigger_status
    #for IceTop
    vem_cal_map = calibration.vem_cal
    station_geo_map = geometry.stationgeo

    bad_dom_list = list()
    if "BadDomsList" in status_frame :
        print("Found a BadDomsList in the frame.  Gonna use it.")
        bad_dom_list = status_frame.Get("BadDomsList")
    else:
        print(status_frame)
        try :
            from icecube.BadDomList import bad_dom_list_static
            bad_dom_list = bad_dom_list_static.IC86_static_bad_dom_list()
        except ImportError :
            print("ERROR : BadDomsList wasn't found in the frame")
            print("and either the BadDomList doesn't exist or")
            print("there's no static_bad_dom_list.")
            sys.exit(1)

    high_qe_dom_list = [icetray.OMKey(36,d) for d in range(44,60) \
               if (d != 45 and d!= 47) ]
    high_qe_dom_list += [ icetray.OMKey( 79, i ) for i in range(30, 45) \
                 if i not in [ 32, 41, 43 ] ]
    high_qe_dom_list += ( [ icetray.OMKey( 80, i ) for i in range(30, 44) ] )
    for s in range(81,87):
        high_qe_dom_list += [icetray.OMKey(s,d) for d in range(1,61)]
        high_qe_dom_list += [icetray.OMKey(43,55)] 

    result_dict = dict()
    result_dict['dom_geo_map'] = dom_geo_map
    result_dict['dom_cal_map'] = dom_cal_map
    result_dict['dom_status_map'] = dom_status_map
    result_dict['trigger_status_map'] = trigger_status_map
    result_dict['vem_cal_map'] = vem_cal_map
    result_dict['station_geo_map'] = station_geo_map
    result_dict['bad_dom_list'] = bad_dom_list
    result_dict['high_qe_dom_list'] = high_qe_dom_list
    
    return result_dict

