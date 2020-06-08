Flux parametrization
====================

MuonGun works by drawing samples from a parametrization of the atmospheric muon
flux as a function of vertical depth, zenith angle, multiplicity, energy, and
for bundles, the distance of each muon from the shower axis. These variables
(the same ones used in a parametrization by `Becherini et al`_) completely
describe a muon bundle under the following assumptions:

1. The flux of cosmic-ray primaries that reach the atmosphere and their 
   daughter muons is independent of the azimuthal arrival direction.
2. The bundle is azimuthally symmetric around the shower axis.
3. All muons in the bundle are perfectly parallel to each other and to the 
   shower axis.
4. The shower front has no curvature.

The first assumption is violated by deflection in the Earth's magnetic field,
but the same approximation is used whenever dCORSIKA showers are randomized in
azimuth before being fed in to IceTray (i.e. nearly always). The remaining
three approximations are important only if it is possible to measure very
detailed properties of the bundle structure over relatively short (~1 km)
observation baselines.

The parameterizations are made by propagating muons from dCORSIKA simulation
filling their properties into a set of histograms
(resources/scripts/propagate_and_fill.py). The distributions in these
histograms are then fitted with tensor-product B-spline surfaces using
photospline (resources/fitting/fit.sh).

.. _`Becherini et al`: http://dx.doi.org/10.1016/j.astropartphys.2005.10.005
