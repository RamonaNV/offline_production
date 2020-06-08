#!/usr/bin/env python

from I3Tray import I3Units
from icecube import icetray
from icecube import dataclasses
from icecube import phys_services
from icecube import cmc

r = phys_services.I3GSLRandomService(23865910)
c = cmc.I3CascadeMCService(r)

nstochastics_l = list()
for i in range(10000):
    particle = dataclasses.I3Particle()

    x = 0.*I3Units.m
    y = 0.*I3Units.m
    z = 0.*I3Units.m
    t = 0.*I3Units.s
    zenith = 0.*I3Units.deg
    azimuth = 270.*I3Units.deg
    e = 40000.*I3Units.GeV

    particle.type = dataclasses.I3Particle.ParticleType.Hadrons
    particle.location_type = dataclasses.I3Particle.LocationType.InIce
    particle.pos = dataclasses.I3Position(x, y, z)
    particle.time = t
    particle.dir = dataclasses.I3Direction(zenith,azimuth)
    particle.energy = e

    daughters = c.Propagate(particle)

    nstochastics_l.append(len(daughters))
    
try:
    import pylab

    pylab.figure()
    pylab.title("Number of Stochastics")
    pylab.xlabel("N")
    pylab.hist(nstochastics_l,histtype="step", log=True)

    pylab.show()

except ImportError:
    print("pylab not installed.  not going to bother plotting anything.")
