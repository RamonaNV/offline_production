#!/usr/bin/env python

"""
read flasher data
needs gcd generated with
gcdserver/resources/BuildGCD.py -r $run -o GCD.i3
"""

from I3Tray import *
from icecube import dataio, dataclasses, payload_parsing
from icecube import DomTools, WaveCalibrator, wavedeform

tray = I3Tray()

tray.Add(dataio.I3Reader, filenamelist=sys.argv[1:])

tray.AddModule("QConverter", WritePFrame=False)

tray.AddSegment(payload_parsing.I3DOMLaunchExtractor,
		MinBiasID = "MinBias",
		FlasherDataID = "Flasher",
		CPUDataID = "BeaconHits",
		SpecialDataID = "SpecialHits",

		# location of scintillators and IceACT
		SpecialDataOMs = [OMKey(0,1),
				  OMKey(12,65),
				  OMKey(12,66),
				  OMKey(62,65),
				  OMKey(62,66)]
		)

def findLeadingEdge(rawAtwdLst, atwdBinSize, ledLight):
    binNum = len(rawAtwdLst)
    initTime = 0
    if (ledLight > 30):
        threshold = 750

        # get bin index of awtd list first element whose pulse < 750
        index = 0
        for i,value in enumerate(rawAtwdLst):
            if (value < threshold):
                index = i
                break

        # Linear interpolation to determine exact time when threshold passed
        initTime = 0
        if (index > 0) and (index < binNum):
            initTime=((rawAtwdLst[index]-threshold)*(index-1.)+(threshold-rawAtwdLst[index-1])*(index))/(rawAtwdLst[index]-rawAtwdLst[index-1])

            if (initTime < 0):
                initTime = 0
            if (initTime > binNum):
                initTime = binNum - 1

        # Now compute the Leading Edge time
        sr = atwdBinSize
        leadEdgeTime = (4 + initTime) * sr
        return leadEdgeTime

    else:

        # compute the minimum pulse from the atwd vector and its index
        minIndex = 0
        minPulse = rawAtwdLst[0]

        for i,value in enumerate(rawAtwdLst):
            if (value < minPulse):
                minPulse = value
                minIndex = i

        # compute leading edge time
        minTime = float(minIndex)
        sr = atwdBinSize
        leadEdgeTime = (4 + initTime) * sr
        return leadEdgeTime

def flasherdecode(frame):
	cal = frame["I3Calibration"]
	status = frame["I3DetectorStatus"]
	domcal = cal.dom_cal
	domstatus = status.dom_status
	event_header = frame["I3EventHeader"]
	run = event_header.run_id
	subrun = event_header.sub_run_id
	runmap = frame["I3FlasherSubrunMap"][subrun]
	event_state = event_header.state

	flashervect = dataclasses.I3FlasherInfoVect()
	if (frame.Has("InIceFlasher") and event_state == 20):
		flashmap = frame["InIceFlasher"]

		for dom,launches in flashmap:
			for launch in launches:

				thisflasher=dataclasses.I3FlasherInfo()

				inspect.getdoc(launch)

				if str(launch.which_atwd) == "ATWDa":
					chip=0
				elif str(launch.which_atwd) == "ATWDb":
					chip=1

				thisflasher.atwd_bin_size = 1./dataclasses.atwd_sampling_rate(chip,domstatus[dom],domcal[dom])

				thisflasher.flashing_om = dom
				thisflasher.atwd_bin_size = 1./dataclasses.atwd_sampling_rate(chip,domstatus[dom],domcal[dom])
				thisflasher.raw_atwd3 = launch.raw_atwd[3]

				if dom in runmap:
					flashInfo = runmap[dom]
					thisflasher.led_brightness = flashInfo.brightness
					thisflasher.mask = flashInfo.mask
					thisflasher.width = flashInfo.window
					thisflasher.rate = flashInfo.rate

					leadingEdge = findLeadingEdge(launch.raw_atwd[3], thisflasher.atwd_bin_size, thisflasher.led_brightness)

					thisflasher.flash_time = launch.time + leadingEdge - 8.3
					flashervect.append(thisflasher)

		if (len(flashervect) > 0):
			frame["FlasherInfo"]=flashervect

tray.AddModule(flasherdecode, Streams=[icetray.I3Frame.DAQ])


def myfilter(frame):
	if (frame.Has("FlasherInfo")):
		flasherinfo = frame.Get("FlasherInfo")
		if (len(flasherinfo) != 1):
			return False
		else:
			return True
	else:
		return False

tray.AddModule(myfilter, Streams=[icetray.I3Frame.DAQ]) # If = lambda f: myfilter(f)

tray.AddModule("I3LCCleaning",
               InIceInput = "InIceRawData",
               InIceOutput = "InIceRawDataClean")

tray.AddModule("I3WaveCalibrator",
               Launches="InIceRawDataClean")

tray.AddModule("I3PMTSaturationFlagger")

tray.AddModule("I3Wavedeform")


def PulseShift(frame):
     if (frame.Has("WavedeformPulses") and frame.Has("FlasherInfo")):
          pulse_map = frame["WavedeformPulses"]
          flashervect = frame.Get("FlasherInfo")
          f = flashervect[0]
          print "F", f.flashing_om.string, f.flashing_om.om, f.flash_time, f.mask, f.led_brightness, f.width, f.rate
          for om,pulse_series in pulse_map:
               for pulses in pulse_series:
                    print "R", om.string, om.om, pulses.time, pulses.charge

          if (frame.Has("CalibrationErrata")):
              errata = frame["CalibrationErrata"]
              for dom,val in errata:
                for tw in val:
                    print "E", dom.string, dom.om, tw.start, tw.stop

          if (frame.Has("SaturationWindows")):
              satwin = frame["SaturationWindows"]
              for dom,val in satwin:
                for tw in val:
                    print "S", dom.string, dom.om, tw.start, tw.stop

          print

tray.AddModule(PulseShift, Streams=[icetray.I3Frame.DAQ])


# tray.AddModule("Dump")

tray.AddModule("TrashCan")
tray.Execute()
tray.Finish()
