#!/usr/bin/env python
##################################################################################
# A script designed to produce events containing only noise triggers. Contains
# both a normal tray segment (with everything except a random module and an
# I3Writer) and an iceprod module for SimProd generation.
#
# Each event produced will initally be a 100 ms long frame of noise. This will
# be triggered and then cut down using CoincidenceAfterProcessing from trigger-sim
# to create normal-sized events.
#
# @author Michael Larson (mjlarson@nbi.ku.dk)
# @date 15 July 2015
##################################################################################
from icecube import icetray

@icetray.traysegment
def ProduceNoiseTriggers(tray, name, gcd_file, nevents=1, run_id=None):
    """
    ProduceNoiseTriggers: Create events containing only noise and no simulated particle interactions. These
    are needed for low-energy DeepCore and PINGU simulation due to a lower detector threshold. There
    are some nuances to this. All events must be long in order to avoid problems at the edge of simulated
    events where expanding around the triggers could include unsimulated regions in time. All frames are
    initially 100 ms long. These are triggered, then cut down by CoincidenceAfterProcessing to be more
    reasonable lengths.

    :param name: Name to add to all of the modules.
    :param gcd_file: The location and name of the GCD file to use
    :param mjd: Modified Julian Date to put into the I3EventHeader. Start of IC86-1 by default.
    :param runID: The run number to put into the I3EventHeader.
    :param nevents: The number of events to simulate.
    """
    tray.AddModule("I3InfiniteSource",name+"_streams",
                   Prefix=gcd_file,
                   Stream=icetray.I3Frame.DAQ)

    from icecube import vuvuzela
    from icecube.icetray import I3Units
    tray.AddSegment(vuvuzela.AddNoise,name+"_vuvuzela",
            InputName = "",
            OutputName = "I3MCPESeriesMap",
            EndTime = 50*I3Units.millisecond,
            StartTime = -50*I3Units.millisecond,
            )

    from icecube.DOMLauncher.domlauncher import DetectorResponse
    tray.AddSegment(DetectorResponse, name+"_detector",
            dom_config= {"BeaconLaunches":False})

    from icecube import trigger_sim, dataio
    tray.AddSegment(trigger_sim.TriggerSim,name+"_triggers",
                    run_id = run_id,
                    gcd_file = dataio.I3File(gcd_file,'r'))

    # Prepare for the CoincidenceAfterProcessing
    from icecube.dataclasses import I3MCTree
    def fakeTree(frame):
        frame["I3MCTree"] = I3MCTree()

    tray.AddModule(fakeTree, name+"_fakeDatTree",
           Streams=[icetray.I3Frame.DAQ,])
    tray.AddModule(NoiseWeight, name+"_addNoiseWeights",
           NEvents = nevents,
           Streams=[icetray.I3Frame.DAQ,])

    tray.AddModule(CheckEventTime, name+"_timecheck")

    # This is normally for long-frame CORSIKA, but its useful here too.
    # Breaks the long frame into multiple smaller one by magic.
    from icecube import trigger_sim
    tray.AddModule("CoincidenceAfterProcessing", name+"breakIntoFrames",
           MinimumSignalNch = 0)

    # Now return to whatever called me
    return


def NoiseWeight(frame, NEvents):
    """
    Write a weightmap to the frame for the noise triggers. REQUIRED for the CoincidenceAfterProcessing module. This
    should only depend on the number of events, since there's no physics or flux models to worry about to first order.
    There is also a "time buffer", which is meant to take into account edge effects of triggers near the frame edges.

    :param NumEvents: Number of events being produced for this file
    """

    from icecube.dataclasses import I3MapStringDouble, I3MCTree
    from icecube.icetray import I3Units
    weightmap = I3MapStringDouble()
    weightmap["simulated_livetime"] = 100*I3Units.millisecond
    weightmap["time_buffer"] = 2*20*I3Units.microsecond
    weightmap["nevents"] = NEvents
    weightmap["weight"] = 1.0/(weightmap["nevents"]*(weightmap["simulated_livetime"]-weightmap["time_buffer"]))
    frame["noise_weight"] = weightmap
    return


def CheckEventTime(frame, TriggerHierarchy="I3TriggerHierarchy", MCPESeriesMap="I3MCPESeriesMap"):
    """
    If the trigger exists within 20 microseconds of the first or last hit, then we need
    to remove that trigger. Otherwise events will include a non-simulated region in the
    event window and things will just be wrong.

    :param TriggerHierarchy: Name of the I3TriggerHierarchy in the frame.
    :param MCPESeriesMap: Name of the I3MCPESeriesMap to check for earliest/latest hits.
    """

    from icecube.dataclasses import I3EventHeader, I3TriggerHierarchy
    eventtriggers = frame[TriggerHierarchy]
    eventmcpes = frame[MCPESeriesMap]

    # We'll need to clone the trigger hierarchy to make live easier.
    newtriggers = I3TriggerHierarchy(eventtriggers)

    # Get the first and last time for the MCPEs
    hittimes = [hit.time for hit in hitseries for hitseries in eventmcpes.values()]
    earliest = min(hittimes)
    latest = max(hittimes)

    # Loop over the trigger times and remove any triggers within 20 microseconds of
    # the earliest/last hit.
    for trigger in eventtriggers:
        if trigger.time - earliest < 20*I3Units.microsecond or\
            latest - trigger.time < 20*I3Units.microsecond:
            newtriggers.erase(trigger)


    del frame[TriggerHierarchy]
    frame[TriggerHierarchy] = newtriggers
    return


