#!/usr/bin/env python

from icecube import icetray

@icetray.traysegment
def Calibration(tray, name,
    BadDOMsHLC = None,
    BadDOMsSLC = None):

    from icecube import icetray, dataclasses, dataio, phys_services
    from icecube import WaveCalibrator
    from icecube import wavedeform
    from icecube import DomTools

    def isOldListSubsetOfNewList(dynamicList, staticList):
        retval=True

        for key in staticList:
            if key not in dynamicList:
                icetray.logging.log_warn("key {} is not in dynamic (in-frame) BadDOMList".format(key))
                retval=False

        return retval

    def AddBadDOMLists(frame, listHLC, listSLC):
        if "BadDomsList" in frame:
            dynamicList = set(frame["BadDomsList"])

            identicalHLC = isOldListSubsetOfNewList(dynamicList, set(listHLC))
            if not identicalHLC:
                icetray.logging.log_warn("*** REPLACING EXISTING BAD DOM LIST in D-frame (HLC)! ***")
                del frame["BadDomsList"]
                newList = set(listHLC).union(dynamicList)
                icetray.logging.log_warn(
                    "  -> new DOM entries in updated list are: {}".format(newList.difference(dynamicList)))
                frame["BadDomsList"] = dataclasses.I3VectorOMKey(newList)
        else:
            frame["BadDomsList"] = dataclasses.I3VectorOMKey(listHLC)

        if "BadDomsListSLC" in frame:
            dynamicList = set(frame["BadDomsListSLC"])

            identicalSLC = isOldListSubsetOfNewList(dynamicList, set(listSLC))
            if not identicalSLC:
                icetray.logging.log_warn("*** REPLACING EXISTING BAD DOM LIST in D-frame (SLC)! ***")
                del frame["BadDomsListSLC"]

                newList = set(listSLC).union(dynamicList)
                icetray.logging.log_warn(
                    "  -> new DOM entries in updated list are: {}".format(newList.difference(dynamicList)))
                frame["BadDomsListSLC"] = dataclasses.I3VectorOMKey(newList)
        else:
            frame["BadDomsListSLC"] = dataclasses.I3VectorOMKey(listSLC)

    if BadDOMsHLC is not None or BadDOMsSLC is not None:
        tray.AddModule(AddBadDOMLists, name+"_AddBadDOMLists",
            listHLC = BadDOMsHLC,
            listSLC = BadDOMsSLC,
            Streams = [icetray.I3Frame.DetectorStatus])

        # online cleaning
        tray.AddModule("I3DOMLaunchCleaning", name+'_baddomclean',
                   CleanedKeys=BadDOMsSLC,
                   InIceOutput='CleanInIceRawData',
                   IceTopOutput='CleanIceTopRawData'
                   )
    else:
        # online cleaning
        tray.AddModule("I3DOMLaunchCleaning", name+'_baddomclean',
                   CleanedKeys=[],
                   InIceOutput='CleanInIceRawData',
                   IceTopOutput='CleanIceTopRawData'
                   )

    # offline cleaning
    tray.AddModule( 'I3DOMLaunchCleaning', name+'_OfflineLaunchCleaning',
        InIceInput = 'CleanInIceRawData',
        IceTopInput = 'CleanIceTopRawData',
        InIceOutput = 'OfflineCleanInIceRawData',
        IceTopOutput = 'OfflineCleanIceTopRawData',
        FirstLaunchCleaning = False,
        CleanedKeysList = 'BadDomsListSLC'
        )

    tray.AddModule('I3WaveCalibrator', name+'_wavecal',
        Launches='OfflineCleanInIceRawData',
        Waveforms='CalibratedWaveforms',
        Errata='CalibrationErrata',
        WaveformRange='CalibratedWaveformRange',
        )
    tray.AddModule('I3Wavedeform', name+'_wavedeform',
        Waveforms='CalibratedWaveforms',
        WaveformTimeRange='CalibratedWaveformRange',
        Output='OfflinePulses',
        )

    tray.AddModule('I3PMTSaturationFlagger', name+'_flag_zorchers',
        Waveforms="CalibratedWaveforms",
        Output="SaturationTimes")

    # I like having frame objects in there even if they are empty for some frames
    def createEmptyTimeWindows(frame, WindowNames=[]):
        for name in WindowNames:
            if name in frame: continue
            frame[name] = dataclasses.I3TimeWindowSeriesMap()
    tray.AddModule(createEmptyTimeWindows, name+'_createEmptyTimeWindows',
        WindowNames = ["CalibrationErrata", "SaturationTimes"],
        Streams=[icetray.I3Frame.DAQ])

    # convert from map<OMKey, I3TimeWindow> to Vector<OMKey> 
    # (i.e. exclude whole DOMs, not only time ranges)
    def convertToVectorOMKey(frame, inputName, outputName):
        if inputName not in frame: return

        keys = []
        origErrata = frame[inputName]
        for key, window in origErrata:
            keys.append(key)
        newObject = dataclasses.I3VectorOMKey(keys)
        frame[outputName] = newObject
    tray.AddModule(convertToVectorOMKey, name+'_convertToVectorOMKey',
        Streams=[icetray.I3Frame.DAQ],
        inputName="SaturationTimes",
        outputName="SaturatedDOMs")

