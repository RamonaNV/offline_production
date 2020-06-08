#!/usr/bin/env python

from icecube import icetray, dataclasses, dataio
from icecube import phys_services, sim_services, simclasses, MuonGun
from I3Tray import I3Tray

tray = I3Tray()

tray.AddModule('I3InfiniteSource', 'driver')
tray.AddService('I3GSLRandomServiceFactory', 'rng', Seed=1337)

surface  = MuonGun.Cylinder(1600, 800)
flux     = MuonGun.BMSSFlux()
flux.min_multiplicity = 1
flux.max_multiplicity = 1
energies = MuonGun.BMSSEnergyDistribution()
radii    = MuonGun.BMSSRadialDistribution()

model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
flux = model.flux
flux.min_multiplicity = 1
flux.max_multiplicity = 1
radii = model.radius
energies = model.energy

generator = 10000*MuonGun.StaticSurfaceInjector(surface, flux, MuonGun.OffsetPowerLaw(2, 500., energies.min, energies.max), radii)

tray.AddModule('I3MuonGun::GeneratorModule', 'GenerateCosmicRayMuons', Generator=generator)
tray.AddModule('I3MuonGun::WeightCalculatorModule', 'weight',
    Model=MuonGun.BundleModel(flux, radii, energies),
    Generator=generator)

weights = []
def check_weight(frame):
    frame['MCMuon'] = frame['I3MCTree'].primaries[0]
    m = len(frame['I3MCTree'].children(frame['MCMuon']))
    weight = frame['weight'].value
    
    # print weight, m
    weights.append(weight)
    assert weight > 0
tray.Add(check_weight, Streams=[icetray.I3Frame.DAQ])


tray.Execute()

try:
    from numpy.testing import assert_approx_equal
    # Because the energies are drawn from a biased distribution, the weight sum
    # converges to the total rate rather slowly. Ask only for agreement within 
    # 10% to catch major mistakes.
    assert_approx_equal(sum(weights), generator.total_rate, 2, err_msg="weights should sum to total rate")
except ImportError:
    pass
