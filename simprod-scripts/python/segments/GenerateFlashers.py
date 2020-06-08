from os.path import expandvars
from icecube import icetray, dataclasses
from I3Tray import I3Units

def GenerateFlashers(frame,
    FlashString=63,
    FlashDOM=20,
    FlashBrightness=127,
    FlashWidth=124,
    FlashTime=0,
    FlashMask=0b000000000001,
    FlasherInfoName="I3FlasherInfo"):

    outputVect = dataclasses.I3FlasherInfoVect()
    flasherInfo = dataclasses.I3FlasherInfo()
    flasherInfo.flashing_om = icetray.OMKey(FlashString,FlashDOM)
    flasherInfo.led_brightness = FlashBrightness
    flasherInfo.mask = FlashMask
    flasherInfo.width = FlashWidth
    flasherInfo.flash_time = FlashTime
    
    flasherInfo.rate = 0
    flasherInfo.atwd_bin_size = 0
    
    outputVect.append(flasherInfo)
    frame[FlasherInfoName] = outputVect



