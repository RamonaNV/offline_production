import os

import numpy as np

import I3Tray
import icecube.icetray
import icecube.dataclasses
import icecube.clsim 
import icecube.simclasses

source_directory = os.environ['I3_SRC']
file = open(source_directory + "/ppc/resources/ice/dx.dat")
contents = []
for line in file:
    contents += [line.split()]
DOMs = []
for key in contents:
    # i contains string , OM_number , orientation of cable , and error on orientation
    DOMs.append( [icecube.icetray.OMKey( int(key[0]) , int(key[1]) ) , float(key[2])] ) 

class AddCylinder(icecube.icetray.I3Module):

    def __init__(self,context):
        icecube.icetray.I3Module.__init__(self,context)
        self.AddParameter( "TopCylinder" , "Position of top of cylinder" , icecube.dataclasses.I3Position( 0 , 0 , 500 ))
        self.AddParameter( "BottomCylinder" , "Position of the bottom of cylinder" , icecube.dataclasses.I3Position( 0 , 0 , -500 ))
        self.AddParameter( "Radius" , "Radius of Cylinder" , 10)
        self.AddParameter( "CableMapName" , "Frame key of cable map" , "CableMap")
        
    def Configure(self):
        self.top = self.GetParameter('TopCylinder')
        self.bottom = self.GetParameter('BottomCylinder')
        self.radius = self.GetParameter('Radius')
        self.cable_map_name = self.GetParameter('CableMapName')

    def Geometry(self,frame):
        frame['Cable'] = icecube.simclasses.I3ExtraGeometryItemCylinder(self.top , self.bottom , self.radius)
        self.PushFrame(frame)
            
class AddCylinders(icecube.icetray.I3Module):
    '''
    For each DOM in the geometry a cylinder is added right next to each DOM.
    '''
    def __init__(self,context):
        icecube.icetray.I3Module.__init__(self,context)
        self.AddParameter( "CylinderLength", "Length of the cable" , 1.0 )
        self.AddParameter( "CylinderRadius" , "Radius of the cable" , 0.023 ) #The radius of the cables is 23 mm
        self.AddParameter( "CableMapName" , "Frame key of cable map" , "CableMap")
        self.AddParameter( "DOMRadius" , "Radius of the DOM" , 0.1651*I3Tray.I3Units.m)
        
    def Configure(self):
        self.height = self.GetParameter("CylinderLength")
        self.radius = self.GetParameter("CylinderRadius")
        self.cable_map_name = self.GetParameter('CableMapName')
        self.dom_radius = self.GetParameter("DOMRadius")
       
    def Geometry(self,frame):
        geometry = frame["I3Geometry"]
        cable_map = icecube.simclasses.I3CylinderMap()
        for i in DOMs:
            om_key = i[0]
            orientation = i[1] * (np.pi/180.0)
            position_x = geometry.omgeo[om_key].position.x + (self.dom_radius + self.radius) * np.cos( orientation )
            position_y = geometry.omgeo[om_key].position.y + (self.dom_radius + self.radius) * np.sin( orientation )
            position_z = geometry.omgeo[om_key].position.z
            top = icecube.dataclasses.I3Position( position_x , position_y , position_z + self.height/2.0)
            bottom = icecube.dataclasses.I3Position( position_x , position_y , position_z - self.height/2.0)
            cable_map[om_key] = icecube.simclasses.I3ExtraGeometryItemCylinder(top, bottom, self.radius)

        frame[self.cable_map_name] = cable_map
        self.PushFrame(frame)

class AverageShadowFraction(icecube.icetray.I3Module):
    def __init__(self, context):
        icecube.icetray.I3Module.__init__(self, context)

        self.AddParameter("PhotonMapName", "Map name of unshadowed photons", "")
        self.AddParameter("ShadowedPhotonMapName", "Map name of shadowed photons", "")

        self.measured_shadow_fractions = list()
        
    def Configure(self):
        self.shadowed_photons = self.GetParameter("ShadowedPhotonMapName")
        self.photons = self.GetParameter("PhotonMapName")
        
    def DAQ(self, frame):
        '''
        Calculate the ratio of shadowed to unshadowed photons.
        '''
        # the first element should be the number of shadowed photons
        # the second element should be the number of unshadowed photons.
        photon_counts = list()
        for frame_key in [self.shadowed_photons, self.photons]:
            photons = frame[frame_key]
            photon_counts.append(sum([len(i) for i in photons.values()]))

        shadow_fraction = float(photon_counts[0])/float(photon_counts[1])
        self.measured_shadow_fractions.append(shadow_fraction)
        self.PushFrame(frame)

    def Finish(self):
        average_shadow_fraction = sum(self.measured_shadow_fractions)/len(self.measured_shadow_fractions)
        print("Average Shadow Fraction = %.2f" % average_shadow_fraction)
        
