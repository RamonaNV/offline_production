Introduction
============

Detector simulation proceeds in 4 steps:

#. Event generation
#. Particle propagation
#. Photon propagation
#. Detector simulation

:cpp:class:`I3PropagatorModule` implements the second step, turning the
:cpp:class:`I3MCTree` produced by event generators (e.g. CORSIKA, MuonGun,
NuGen) into a set of final, Cherenkov-emitting states (muon tracks with finite
length and electromagnetic cascades). Each particle may produce one or more
daughter particles as it propagates (e.g. energy losses of muons, or muons
produced in hadronic showers), which in turn may produce daughter particles of
their own.

The final :cpp:class:`I3MCTree` may be 10 to 100 times larger than the input
tree. Because the transformation from initial to final tree is a function only
of the initial state, the transformation algorithm, and a sequence of random
numbers, however, it can be deleted after photon propagation and reproduced if
need be. :cpp:class:`I3PropagatorModule` stores the state of the configured
random number generator and will recreate it if necessary.

Details
=======

:cpp:class:`I3PropagatorModule` takes a
:cpp:class:`std::map\<I3Particle::ParticleType, I3PropagatorServicePtr>` that
defines the propagator to use for each particle type. It then iterates over the
input :cpp:class:`I3MCTree`, passing each particle in the tree to the
corresponding propagator, if defined. The propagator may modify the particle
(e.g. by setting its length), and also return secondary particles. These are
attached as daughters of the propagated particle in the output tree. If the
secondaries are also propagatable, they are added to the propagation stack for
later processing. This allows for unlimited chains of decays, for example a tau
lepton that decays to a muon which produces photonuclear losses which then also
produce secondary muons.



