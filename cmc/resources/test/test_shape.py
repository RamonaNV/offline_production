#!/usr/bin/env python

from icecube import cmc, phys_services, dataclasses
import math

rng = phys_services.I3GSLRandomService(0)

mc = cmc.I3CascadeMCService(rng)

p = dataclasses.I3Particle()
p.dir = dataclasses.I3Direction(0,0)
p.pos = dataclasses.I3Position(0,0,0)
p.time = 0
p.energy = 1e9
p.type = p.EMinus
p.shape = p.Cascade

radiation_length = 0.358 / 0.9216

p.type = p.EMinus
daughters = mc.Propagate(p)
for d in daughters:
	try:
		assert d.type == d.EMinus, "all daughters are EMinus"
		assert abs(d.length-3*radiation_length) < 1e-7, "segment is 3 radiation lengths"
		assert d.shape == d.CascadeSegment, "shape is cascade segment"
	except AssertionError:
		print(d)
		raise

p.type = p.Hadrons
for d in daughters:
	try:
		if abs(d.type) == d.MuMinus:
			assert math.isnan(d.length)
		else:
			assert d.type == d.EMinus, "all cascade daughters are EMinus"
			assert abs(d.length-3*radiation_length) < 1e-7, "segment is 3 radiation lengths"
			assert d.shape == d.CascadeSegment, "shape is cascade segment"
	except AssertionError:
		print(d)
		raise

