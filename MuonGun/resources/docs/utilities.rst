Utilities: recovering the energies of muons at depth
====================================================

Muons are ranged particles, so both their position and energy are a function of
time. The position is easy enough to represent in an :cpp:class:`I3Particle`;
by dint of special relativity, the position be approximated by speed-of-light
displacement along the muon's initial direction from an arbitrary reference
point. Information about a muon's true energy at any given moment, however, is
scattered in multiple locations in each frame. These are:

* The initial energy at the reference point, stored in the :cpp:class:`I3Particle`
  that represents the muon itself
* The energy the muon has when it enters and exits the MMC simulation volume
  (typically a 1600 x 800 meter upright cylinder), stored in the :cpp:class:`I3MMCTrack`
  corresponding to the muon.
* The sizes of the stochastic energy losses inside the MMC simulation volume. These
  are stored as :cpp:class:`I3Particle` attached as daughters of the :cpp:class:`I3Particle`
  represending the muon in the :cpp:class:`I3MCTree`

.. py:currentmodule:: icecube.MuonGun

While it is in principle straightforward to recover the energy of a simulated
muon at any point along its path inside the MMC simulation volume from these
scattered data structures, it is also extremely tedious. MuonGun includes
a utility class, :py:class:`Track`, to automate this task.

.. py:class:: Track
   :noindex:

   A subclass of I3Particle that includes the particle's energy losses

   .. classmethod:: harvest(frame)
      :noindex:

      Assemble a collection of :py:class:`Track` from the
      :cpp:class:`MMCTrackList` and :cpp:class:`I3MCTree` in *frame*

   .. method:: get_energy(displacement)
      :noindex:

      Get the energy of the muon *displacement* meters from its starting position

This class can be used, for example, to find the total energy losses
(stochastic and continuous) of all muons within some volume::

    from icecube import MuonGun, simclasses
 
    # A surface approximating the actual detector (make it smaller if you only care e.g. about DeepCore)
    surface = MuonGun.Cylinder(1000,500)

    edep = 0
    for track in MuonGun.Track.harvest(frame['I3MCTree'], frame['MMCTrackList']):
        # Find distance to entrance and exit from sampling volume 
        intersections = surface.intersection(track.pos, track.dir)
        # Get the corresponding energies
        e0, e1 = track.get_energy(intersections.first), track.get_energy(intersections.second)
        # Accumulate
        edep +=  (e0-e1)

There is also a convenience function, :py:func:`muons_at_surface`, that uses
:py:class:`Track` to "shift" all the muons in a frame forward so that their
positions, times, and energies correspond to their intersections with a given
surface.

.. py:function:: muons_at_surface(frame, surface)
   :noindex:
 
   :param frame: an I3Frame
   :param surface: a :py:class:`Surface`
   :returns: a list of I3Particles with positions, times, and energies that
             correspond to their intersection with *surface*. Muons that range
             out before reaching *surface* are not included.

This can be used to quickly obtain the true multiplicity of a muon bundle when
it enters the detector::

    surface = MuonGun.Cylinder(1000,500)
    multiplicity = len(MuonGun.muons_at_surface(frame, surface))
