#!/usr/bin/env python

from __future__ import print_function
from icecube import dataclasses, phys_services, PROPOSAL
from math import acos, pi

p = dataclasses.I3Particle(dataclasses.I3Position(0,0,-850), dataclasses.I3Direction(0,0,1), 0)
p.location_type = p.InIce
p.type = p.TauMinus
p.energy = 1e5

prop = PROPOSAL.I3PropagatorServicePROPOSAL()
prop.SetRandomNumberGenerator(phys_services.I3GSLRandomService(0))

products = []
for _ in range(20):
    p.length = 0
    daughters = prop.Propagate(p)
    products.append(daughters)

try:
    assert any([p.MuMinus in [pp.type for pp in daughters] for daughters in products]), "taus decay to muons even outside active volume"
except AssertionError:
    print([[pp.type for pp in daughters] for daughters in products])
    raise

for daughters in products:
    try:
        mu = [pp for pp in daughters if pp.type == pp.MuMinus][0]
    except IndexError:
        continue
    psi = 180*(acos(mu.dir*p.dir))/pi
    assert psi < 10, "muon is roughly parallel to parent ({:.1f} < 10 degrees)".format(psi)