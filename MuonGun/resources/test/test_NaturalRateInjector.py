#!/usr/bin/env python

from icecube import icetray, dataclasses, dataio
from icecube import phys_services, sim_services, simclasses, MuonGun
from I3Tray import I3Tray

tray = I3Tray()

tray.AddModule('I3InfiniteSource', 'driver')
tray.AddService('I3GSLRandomServiceFactory', 'rng', Seed=1337)

surface  = MuonGun.Cylinder(1600, 800)
model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
flux = model.flux
radii = model.radius
energies = model.energy

nevents = 10
generator = nevents*MuonGun.NaturalRateInjector(surface, flux, energies)

tray.AddModule('I3MuonGun::GeneratorModule', 'generator', Generator=generator)
tray.AddModule('I3MuonGun::WeightCalculatorModule', 'weight',
    Model=MuonGun.BundleModel(flux, radii, energies),
    Generator=generator)

weights = []
def check_weight(frame):
    frame['MCMuon'] = frame['I3MCTree'].primaries[0]
    m = len(frame['I3MCTree'].children(frame['MCMuon']))
    weight = frame['weight'].value
    
    # print weight, m, generator.total_rate, generator.total_rate/nevents/weight
    weights.append(weight)
    assert weight > 0
tray.Add(check_weight, Streams=[icetray.I3Frame.DAQ])

tray.AddModule('TrashCan', 'YesWeCan')
tray.Execute()
tray.Finish()

try:
    from numpy.testing import assert_approx_equal
    # These should match by construction, modulo round-off
    assert_approx_equal(sum(weights), generator.total_rate, 12, err_msg="weights should sum to total rate")
except ImportError:
    pass