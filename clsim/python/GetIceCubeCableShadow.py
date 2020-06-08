
from icecube.icetray import OMKey
from icecube.simclasses import I3MapModuleKeyI3ExtraGeometryItemCylinder, I3ExtraGeometryItemCylinder
from icecube.dataclasses import I3Position, ModuleKey

from I3Tray import I3Units

import numpy as np
from os.path import expandvars

from_cable_shadow = expandvars("$I3_BUILD/ice-models/resources/models/cable_position/orientation.cable_shadow.txt")
from_led7 = expandvars("$I3_BUILD/ice-models/resources/models/cable_position/orientation.led7.txt")

def GetIceCubeCableShadow(CableAngles=from_led7,
    DOMRadius=165.1*I3Units.mm, CableRadius=23*I3Units.mm, CableLength=1*I3Units.m):
    """
    Get a cylinder representing the position of the cable at each DOM
    
    :param CableAngles: text file containing string, om, angle (degrees), angle error (degrees)
    :param DOMRadius: radius of the DOM sphere
    :param CableRadius: radius of the cable
    :param CableLength: length of the cable segment at each DOM
    :returns: a map of I3ExtraGeometryItem representing the local position of
        the cable *in DOM-centered coordinates*
    """
    # assume the cable runs along the surface of the DOM
    radius = DOMRadius + CableRadius
    shadows = I3MapModuleKeyI3ExtraGeometryItemCylinder()
    for string, om, angle, _ in np.loadtxt(CableAngles, dtype=[('string',int),('om',int),('angle',float),('angle_err',float)]):
        pos = I3Position(radius*np.cos(np.radians(angle)), radius*np.sin(np.radians(angle)), 0)
        shadows[ModuleKey(int(string),int(om))] = I3ExtraGeometryItemCylinder(pos + I3Position(0,0,CableLength/2.), pos + I3Position(0,0,-CableLength/2.), CableRadius)

    return shadows
