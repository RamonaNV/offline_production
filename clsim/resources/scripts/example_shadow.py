#!/usr/bin/env python 

import sys
import os
import os.path
import logging

import numpy

import icecube.icetray
import icecube.dataio
import icecube.clsim
import icecube.clsim.shadow
import icecube.clsim.shadow.cylinder_utils
import icecube.dataclasses
import icecube.phys_services
import icecube.simclasses

from I3Tray import I3Tray
from I3Tray import I3Units

randomService = icecube.phys_services.I3GSLRandomService(seed = 10)
gcd_file = os.getenv('I3_TESTDATA') + '/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'
cable_orientations = os.path.expandvars("$I3_BUILD/ice-models/resources/models/cable_position/orientation.led7.txt")

class InjectCables(icecube.icetray.I3Module):
    '''
    This module injects cylinders as cables reading the measured orientations
    stored in a text file in the ice-models project.
    '''
    def __init__(self, context):
        icecube.icetray.I3Module.__init__(self, context)
        self.AddParameter("CableMapName", "Frame key of the cable map",'CableMap')

    def Configure(self):
        self.cable_map_name = self.GetParameter("CableMapName")
        
    def Geometry(self, frame):
        dom_radius=165.1*I3Units.mm
        cable_radius=23*I3Units.mm
        cable_length=10*I3Units.m

        radius = dom_radius + cable_radius
        cylinder_map = icecube.simclasses.I3CylinderMap()

        geo = frame['I3Geometry'].omgeo
        
        for string, om, angle, _ in numpy.loadtxt(cable_orientations,
                                                  dtype=[('string',int),('om',int),('angle',float),('angle_err',float)]):
            dx = radius*numpy.cos(numpy.radians(angle))
            dy = radius*numpy.sin(numpy.radians(angle))
            
            top_relative = icecube.dataclasses.I3Position(dx, dy, cable_length/2.)
            bottom_relative = icecube.dataclasses.I3Position(dx, dy, -cable_length/2.)
            
            omkey = icecube.icetray.OMKey(int(string),int(om))
            dom_position = geo[omkey].position

            # want positions in lab/detector coordinates
            cylinder_top = dom_position + top_relative
            cylinder_bottom = dom_position + bottom_relative
            cylinder_map[omkey] = icecube.simclasses.I3ExtraGeometryItemCylinder(cylinder_top, cylinder_bottom, cable_radius)
                
        frame[self.cable_map_name] = cylinder_map
        self.PushFrame(frame)        
        
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


tray = I3Tray()
cable_map_name = 'CableMap'
tray.Add("I3InfiniteSource" , 
         Prefix=os.path.expandvars(gcd_file),
         Stream=icecube.icetray.I3Frame.DAQ)

# AddCylinders puts a cylinder next to each DOM using
# the same orientation for each cylinder.
# It's recommended to use the measured orientations,
# which is what InjectCables does.
#tray.Add(icecube.clsim.shadow.cylinder_utils.AddCylinders,
#         CableMapName = cable_map_name,
#         CylinderLength = 17.0,
#         CylinderRadius = 0.023)
tray.Add(InjectCables, CableMapName=cable_map_name)

tray.Add("Dump")

tray.Add(GenerateEvent,
         Type = icecube.dataclasses.I3Particle.ParticleType.MuMinus,
         NEvents = 10,
         XCoord = -256.14,
         YCoord = -521.08,
         ZCoord = 496.03,
         PrimaryDirection = icecube.dataclasses.I3Direction(0 , 0 ,-1),
         DaughterDirection = icecube.dataclasses.I3Direction(0 , 0 , -1),
         I3RandomService = randomService,
         Energy = 10.0*I3Units.TeV )

photonSeriesName = "Photons"
usegpus = any([device.gpu for device in icecube.clsim.I3CLSimOpenCLDevice.GetAllDevices()])    
tray.Add(icecube.clsim.I3CLSimMakePhotons,
         UseGPUs = usegpus,
         UseOnlyDeviceNumber=0,
         UseCPUs = not usegpus,                    
         PhotonSeriesName = photonSeriesName,
         RandomService = randomService,
         IceModelLocation = os.path.expandvars("$I3_BUILD/ice-models/resources/models/spice_lea"),
         GCDFile = gcd_file)

tray.Add("I3ShadowedPhotonRemoverModule",
         InputPhotonSeriesMapName = photonSeriesName,
         OutputPhotonSeriesMapName = photonSeriesName+'Shadowed',
         CableMapName = cable_map_name,
         Distance = 20.0)

tray.AddModule(icecube.clsim.shadow.cylinder_utils.AverageShadowFraction,
               PhotonMapName=photonSeriesName,
               ShadowedPhotonMapName=photonSeriesName+'Shadowed',
)

tray.Add("I3Writer", filename='shadowed_photons_removed.i3')

tray.Execute()
