#!/usr/bin/env python 

import sys
import os
import os.path
import unittest

from I3Tray import I3Tray
from I3Tray import I3Units

import icecube.icetray
import icecube.dataio
import icecube.clsim
import icecube.clsim.shadow.cylinder_utils
import icecube.dataclasses
import icecube.phys_services
import icecube.simclasses

photonSeriesMapName = "Photons"
randomService = icecube.phys_services.I3GSLRandomService(seed = 1)
gcd_file = os.getenv('I3_TESTDATA') + '/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'

class GenerateEvent(icecube.icetray.I3Module):
    def __init__(self, context):
        icecube.icetray.I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("Type", "", icecube.dataclasses.I3Particle.ParticleType.EMinus)
        self.AddParameter("Energy", "", 10.*I3Units.TeV)
        self.AddParameter("NEvents", "", 1)
        self.AddParameter("XCoord", "", 0.)
        self.AddParameter("YCoord", "", 0.)
        self.AddParameter("ZCoord", "", 0.)
        self.AddParameter("PrimaryDirection","",icecube.dataclasses.I3Direction(0.,0.,-1.))
        self.AddParameter("DaughterDirection","",icecube.dataclasses.I3Direction(0.,0.,-1.))
        self.AddOutBox("OutBox")

    def Configure(self):
        self.rs = self.GetParameter("I3RandomService")
        self.particleType = self.GetParameter("Type")
        self.energy = self.GetParameter("Energy")
        self.nEvents = self.GetParameter("NEvents")
        self.xCoord = self.GetParameter("XCoord")
        self.yCoord = self.GetParameter("YCoord")
        self.zCoord = self.GetParameter("ZCoord")
        self.primary_direction = self.GetParameter("PrimaryDirection")
        self.daughter_direction = self.GetParameter("DaughterDirection")

        self.eventCounter = 0

    def DAQ(self, frame):
        daughter = icecube.dataclasses.I3Particle()
        daughter.type = self.particleType
        daughter.energy = self.energy
        daughter.pos = icecube.dataclasses.I3Position(self.xCoord,self.yCoord,self.zCoord)
        daughter.dir = self.daughter_direction
        daughter.time = 0.
        daughter.location_type = icecube.dataclasses.I3Particle.LocationType.InIce

        primary = icecube.dataclasses.I3Particle()
        primary.type = icecube.dataclasses.I3Particle.ParticleType.NuE
        primary.energy = self.energy
        primary.pos = icecube.dataclasses.I3Position(self.xCoord,self.yCoord,self.zCoord)
        primary.dir = self.primary_direction
        primary.time = 0.
        primary.location_type = icecube.dataclasses.I3Particle.LocationType.Anywhere

        mctree = icecube.dataclasses.I3MCTree()
        mctree.add_primary(primary)
        mctree.append_child(primary,daughter)

        frame["I3MCTree"] = mctree

        self.PushFrame(frame)

        self.eventCounter += 1
        if self.eventCounter==self.nEvents:
            self.RequestSuspension()

class SanityCheck(unittest.TestCase):
    photons = photonSeriesMapName
    shadowed_photons = photonSeriesMapName+"Shadowed"
    
    def testKeys(self):
        self.assertTrue(self.shadowed_photons in self.frame, "The shadowed_photons actually shows up in the frame.")
        
    def testShadowRemoval(self):
            
        non_shadowed = self.frame[self.photons]
        count_photons = 0
        
        for i in non_shadowed.values():
            count_photons += len(i)
                
        shadowed = self.frame[self.shadowed_photons]
        count_shadowed = 0
        for i in shadowed.values():
            count_shadowed+= len(i)
                
        icecube.icetray.logging.log_info("photons before shadow: {} after: {}".format(count_photons, count_shadowed),
                                         unit="testCableShadow")

        self.assertGreater(count_photons, 0, "some photons made it to DOMs")
        self.assertGreater(count_photons, count_shadowed, "some photons removed")
        self.assertAlmostEqual(count_shadowed/float(count_photons), 0.96, delta=0.03, msg="approximately 10% photons removed")
            
if __name__ == "__main__":
    icecube.icetray.logging.set_level_for_unit("testShadowPhotonRemover", "INFO")

    tray = I3Tray()

    cable_map_name = 'CableMap'
    
    tray.AddModule("I3InfiniteSource" ,
                   Prefix=os.path.expandvars(gcd_file),
                   Stream = icecube.icetray.I3Frame.DAQ)

    tray.Add(icecube.clsim.shadow.cylinder_utils.AddCylinders,
             CableMapName = cable_map_name,
             CylinderLength = 17.0,
             CylinderRadius = 0.023)

    tray.AddModule("I3MCEventHeaderGenerator",
                   Year = 2009,
                   DAQTime=158100000000000000,
                   RunNumber=1,
                   EventID=1,
                   IncrementEventID=True)

    tray.AddModule(GenerateEvent, 
                   Type = icecube.dataclasses.I3Particle.ParticleType.EMinus,
                   NEvents = 1000,
                   XCoord = -256.14,
                   YCoord = -521.08,
                   ZCoord = 496.03,
                   PrimaryDirection = icecube.dataclasses.I3Direction(0 , 0 ,-1),
                   DaughterDirection = icecube.dataclasses.I3Direction(0 , 0 , -1),
                   I3RandomService = randomService,
                   Energy = 1.0*I3Units.TeV )

    # prefer GPUs
    usegpus = any([device.gpu for device in icecube.clsim.I3CLSimOpenCLDevice.GetAllDevices()])    
    tray.AddSegment(icecube.clsim.I3CLSimMakePhotons,
                    UseGPUs = usegpus,
                    UseOnlyDeviceNumber=0,
                    UseCPUs = not usegpus,                    
                    PhotonSeriesName = photonSeriesMapName,
                    MCPESeriesName = None,
                    RandomService = randomService,
                    IceModelLocation = os.path.expandvars("$I3_BUILD/ice-models/resources/models/spice_lea"),
                    CableOrientation = None,
                    GCDFile = gcd_file)
    
    tray.Add("I3ShadowedPhotonRemoverModule",
             InputPhotonSeriesMapName = photonSeriesMapName,
             OutputPhotonSeriesMapName = photonSeriesMapName+'Shadowed',
             CableMapName = cable_map_name,
             Distance = 10*I3Units.m)
                
    tray.AddModule(icecube.icetray.I3TestModuleFactory(SanityCheck),
                   streams=[icecube.icetray.I3Frame.DAQ])

    tray.AddModule(icecube.clsim.shadow.cylinder_utils.AverageShadowFraction,
                   PhotonMapName=photonSeriesMapName,
                   ShadowedPhotonMapName=photonSeriesMapName+'Shadowed',
    )
    
    tray.Execute()
    
