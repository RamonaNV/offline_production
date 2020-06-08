Rationale
=========

In diffuse neutrino searches it is almost always necessary to estimate the
background due to atmospheric muons from simulation. In IceCube this has
typically been done directly, by simulating air showers to ground level with
CORSIKA, propagating the muons in the shower through the firn and ice with
PROPOSAL to a cylindrical sampling surface surrounding the detector, and then
weighting the simulated events to an assumed cosmic-ray flux. Though this
method offers the highest possible precision available from the chosen
simulation software, it suffers from two key inefficiencies. First, since
the simulation starts with cosmic-ray primaries rather than in-ice muons,
one has only loose control over the characteristics of the muon bundles that
actually reach the detector. For example, if one were interested only in
single muons with a few TeV of energy, one would spend quite a lot of time
simulating both showers whose muons never reached the detector those that
result in high-multiplicity bundles. Second, the direct approach makes it
necessary to repeat the entire simulation chain in order to change aspects
of the air shower simulation such as atmospheric profile or hadronic model.

An alternative approach is to de-couple the air shower simulation and muon
propagation from the remainder of the simulation by constructing a
parametrization of the muon flux under the ice and drawing muon bundles from
the parameterized distribution. This allows one to generate specific bundle
configurations and weight them properly, and also to re-weight existing
simulated events to a muon flux associated with different assumptions about
interactions in the atmosphere.

The parametric approach is used heavily by ANTARES in the form of
their `MUPAGE`_ event generator. The work described here is an attempt to
apply the technique, described in a paper by `Becherini et al`_, to IceCube
simulation.

.. _MUPAGE: http://dx.doi.org/10.1016/j.cpc.2008.07.014
.. _`Becherini et al`: http://dx.doi.org/10.1016/j.astropartphys.2005.10.005
