'''
A collection of functions that are useful shortcuts for pulling
things out of GCD files.
'''
from icecube import icetray, dataclasses, dataio

def get_time(gcdfile):
    '''
    Gets the start time of the I3DetectorStatus frame.
    This was used to generate a DrivingTime from a GCD, so
    is likely obsolete now.
    '''
    frame = gcdfile.pop_frame()
    while not frame.Has("I3DetectorStatus"):
        frame = gcdfile.pop_frame()

    return frame.Get("I3DetectorStatus").start_time

def get_omgeo(gcdfile):
    '''
    Gets the I3OMGeoMap from a GCD file.
    '''
    frame = gcdfile.pop_frame()
    while not frame.Has("I3Geometry"):
        frame = gcdfile.pop_frame()

    return frame.Get("I3Geometry").omgeo

def get_domcal(gcdfile):
    '''
    Gets the I3DOMCalibrationMap from a GCD file.
    '''
    frame = gcdfile.pop_frame()
    while not frame.Has("I3Calibration"):
        frame = gcdfile.pop_frame()

    return frame.Get("I3Calibration").dom_cal

def get_domstatus(gcdfile):
    '''
    Gets the I3DOMStatusMap from a GCD file.
    '''
    frame = gcdfile.pop_frame()
    while not frame.Has("I3DetectorStatus"):
        frame = gcdfile.pop_frame()

    return frame.Get("I3DetectorStatus").dom_status

def get_triggerstatus(gcdfile):
    '''
    Gets the I3TriggerStatus from a GCD file.
    '''
    frame = gcdfile.pop_frame()
    while not frame.Has("I3DetectorStatus"):
        frame = gcdfile.pop_frame()

    return frame.Get("I3DetectorStatus").trigger_status

def put_triggerstatus(trigger_status_map,
                      gcdfile,
                      output_gcd_filename):
    '''
    Adds an I3TriggerStatus to a GCD file.
    '''
    newgcd = dataio.I3File(output_gcd_filename, dataio.I3File.Mode.Writing)
    frame = gcdfile.pop_frame()

    while not frame.Has("I3DetectorStatus"):
        newgcd.push(frame)
        frame = gcdfile.pop_frame()
        
    ds = frame.Get("I3DetectorStatus")
    del frame["I3DetectorStatus"]
    ds.trigger_status = trigger_status_map
    frame["I3DetectorStatus"] = ds
    newgcd.push(frame)

    return newgcd
