#!/usr/bin/env python 

import sys
import os
import os.path

import I3Tray
import icecube.icetray
import icecube.clsim
import icecube.dataclasses
import icecube.phys_services

randomService = icecube.phys_services.I3GSLRandomService(seed = 1)
gcd_file = os.getenv('I3_TESTDATA') + '/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'

class TestShadowFraction(icecube.icetray.I3Module):
    def __init__(self, context):
        icecube.icetray.I3Module.__init__(self, context)

        self.AddParameter("PhotonMapName", "Frame key of unshadowed photons", "Photons")
        self.AddParameter("ShadowedPhotonMapName", "Frame key of shadowed photons", "PhotonsAfterShadow")

        self.shadow_fractions = list()

    def Configure(self):
        self.photon_map_name = self.GetParameter("PhotonMapName")
        self.shadowed_photon_map_name = self.GetParameter("ShadowedPhotonMapName")
        
    def DAQ(self, frame):
        self.PushFrame(frame)
        photon_counts = list()
        for frame_key in [self.shadowed_photon_map_name, self.photon_map_name]:
            photons = frame[frame_key]
            photon_counts.append(sum([len(i) for i in photons.values()]))

        shadow_fraction = float(photon_counts[0])/float(photon_counts[1])
        self.shadow_fractions.append(shadow_fraction)

    def Finish(self):
        # histogram and fit.
        try:
            import numpy
            import scipy.optimize

            bin_values, bin_edges = numpy.histogram(self.shadow_fractions)
            bin_width = bin_edges[1] - bin_edges[0]
            bin_centers = [edge + bin_width/2. for edge in bin_edges[:-1]]
            gaussian = lambda x, A, x0, sigma: A*numpy.exp(-0.5*((x-x0)/sigma)**2) # Gaussian
            p0 = [250, 0.9, 0.1] # Initial guess for the parameters
            popt, pcov = scipy.optimize.curve_fit(gaussian, bin_centers, bin_values, p0)
            print("Fit Result = %s " % str(popt))
            print("Covariance = %s " % str(pcov))

            mean_fit = popt[1]
            assert(mean_fit > 0.89 and mean_fit < 0.925)
            
            # Uncomment the following lines to render
            # the histogram and fit result
            # 
            #import pylab            
            #pylab.hist(self.shadow_fractions)
            #X = numpy.arange(bin_centers[0], bin_centers[-1], 0.01)
            #Y = gaussian(X, popt[0], popt[1], popt[2])
            #pylab.plot(X, Y, '-')
            #pylab.show()
        except ImportError:
            average_shadow_fractions = sum(self.shadow_fractions)/float(len(self.shadow_fractions))
            print(average_shadow_fractions)
            assert(average_shadow_fractions > 0.85 and average_shadow_fractions < 0.95)

            
            
class GenerateEvent(icecube.icetray.I3Module):
    def __init__(self, context):
        icecube.icetray.I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("Type", "", icecube.dataclasses.I3Particle.ParticleType.EMinus)
        self.AddParameter("Energy", "", 10.*I3Tray.I3Units.TeV)
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
        daughter.pos = icecube.dataclasses.I3Position(self.xCoord, self.yCoord, self.zCoord)
        daughter.dir = self.daughter_direction
        daughter.time = 0.
        daughter.location_type = icecube.dataclasses.I3Particle.LocationType.InIce

        primary = icecube.dataclasses.I3Particle()
        primary.type = icecube.dataclasses.I3Particle.ParticleType.NuE
        primary.energy = self.energy
        primary.pos = icecube.dataclasses.I3Position(self.xCoord, self.yCoord, self.zCoord)
        primary.dir = self.primary_direction
        primary.time = 0.
        primary.location_type = icecube.dataclasses.I3Particle.LocationType.Anywhere

        mctree = icecube.dataclasses.I3MCTree()
        mctree.add_primary(primary)
        mctree.append_child(primary, daughter)

        frame["I3MCTree"] = mctree

        self.PushFrame(frame)

        self.eventCounter += 1
        if self.eventCounter==self.nEvents:
            self.RequestSuspension()

if __name__ == "__main__":
    icecube.icetray.logging.set_level_for_unit("testCableShadow", "INFO")

    tray = I3Tray.I3Tray()

    tray.AddModule("I3InfiniteSource" ,
                   Prefix=os.path.expandvars(gcd_file),
                   Stream = icecube.icetray.I3Frame.DAQ)

    tray.AddModule("I3MCEventHeaderGenerator",
                   Year=2009,
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
                   Energy = 1.0*I3Tray.I3Units.TeV )

    photonSeriesName = "Photons"
    # prefer GPUs
    usegpus = any([device.gpu for device in icecube.clsim.I3CLSimOpenCLDevice.GetAllDevices()])    
    tray.AddSegment(icecube.clsim.I3CLSimMakePhotons,"MakePhotons",
                    UseGPUs = usegpus,
                    UseOnlyDeviceNumber=0,
                    UseCPUs = not usegpus,                    
                    PhotonSeriesName = photonSeriesName,
                    MCPESeriesName = None,
                    RandomService = randomService,
                    IceModelLocation = os.path.expandvars("$I3_BUILD/ice-models/resources/models/spice_lea"),
                    CableOrientation = None,
                    GCDFile = gcd_file)
    tray.AddSegment(icecube.clsim.I3CLSimMakePhotons,"MakePhotonsWithShadow",
                    UseGPUs = usegpus,
                    UseOnlyDeviceNumber=0,
                    UseCPUs = not usegpus,
                    PhotonSeriesName = photonSeriesName+"AfterShadow",
                    MCPESeriesName = None,
                    RandomService = randomService,
                    IceModelLocation = os.path.expandvars("$I3_BUILD/ice-models/resources/models/spice_lea"),
                    CableOrientation = os.path.expandvars("$I3_BUILD/ice-models/resources/models/cable_position/orientation.led7.txt"),
                    GCDFile = gcd_file)

    tray.AddModule(TestShadowFraction)    

    tray.Execute()
    
