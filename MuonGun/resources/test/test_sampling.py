#!/usr/bin/env python

"""
Test that samples from the joint radius/energy distribution are actually
distributed like the model
"""

from icecube import icetray, MuonGun, phys_services

import numpy
from scipy import integrate
import warnings
import resource

rng = phys_services.I3GSLRandomService(0)

model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
# only look at bundles with at least 10 muons
model.flux.min_multiplicity = 10

def draw_bundle(rng, flux):
    """
    Generate a bundle impact point, zenith angle, and multiplicity proportional
    to the total flux
    """
    max_flux = float(flux(1.5, 1., flux.min_multiplicity))
    while True:
        depth = rng.uniform(1.5, 2.5)
        ct = rng.uniform(0, 1)
        m = flux.min_multiplicity + rng.integer(flux.max_multiplicity-flux.min_multiplicity)
        if rng.uniform(0, max_flux) < float(flux(depth, ct, m)):
            return depth, ct, m

def sample_energy(edist, depth, ct, m, nsamples=10000):
    
    # bins that will have roughly equal contents
    rbins = numpy.array([0, 1, 2, 3, 4, 6, 10, 15, 250])
    powerlaw = MuonGun.OffsetPowerLaw(4, 1e3, edist.min, edist.max)
    ebins = powerlaw.isf(numpy.linspace(0, 1, 21)[::-1])
    
    start = resource.getrusage(resource.RUSAGE_SELF)
    samples = edist.generate(rng, depth, ct, m, nsamples)
    end = resource.getrusage(resource.RUSAGE_SELF)
    dt = end.ru_utime - start.ru_utime
    icetray.logging.log_info("%.1f microsec/sample" % (1e6*dt/nsamples))
    samples = numpy.array([(i.first, i.second) for i in samples])
    
    bins, edges = numpy.histogramdd(samples, bins=(rbins, ebins))
    assert bins.sum() == nsamples
    
    empties = (bins < 10).sum()/float(bins.size)
    if empties > 0.25:
        warnings.warn('%d%% of bins have fewer than 10 entries' % (100*empties))
    
    norm = nsamples/edist.integrate(depth, ct, m, 0, 250, edist.min, edist.max)
    @numpy.vectorize
    def integrate_model(rbin, ebin):
        return edist.integrate(depth, ct, m, edges[0][rbin], edges[0][rbin+1], edges[1][ebin], edges[1][ebin+1])
    
    i, j = numpy.meshgrid(range(len(rbins)-1), range(len(ebins)-1), indexing='ij')
    mu = norm*integrate_model(i.T, j.T)
    chi2 = (bins.T - mu)**2/mu
    return samples, chi2.sum(), bins.size

for i in xrange(10):
    depth, ct, m = draw_bundle(rng, model.flux)
    samples, chi2, ndof = sample_energy(model.energy, depth, ct, m, nsamples=10000)
    icetray.logging.log_info('depth=%.1f kmwe, cos(theta)=%.2f, m=%d, chi2/ndof = %.2f' % (depth, ct, m, chi2/ndof))
    assert chi2/ndof < 1.5, "Samples follow model"
