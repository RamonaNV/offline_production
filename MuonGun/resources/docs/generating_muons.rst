Generating muons
================

Muon bundle generators are instances of the :cpp:class:`Generator` class. There are currently 3 generators implemented:

.. cpp:class:: StaticSurfaceInjector

   Injects muon bundles on a fixed sampling surface with depth and zenith
   angle distributions given by a flux model and an energy distribution
   given by a broken power law. This is more or less identical to `MUPAGE`_.

.. cpp:class:: EnergyDependentSurfaceInjector

   Similar to :cpp:class:`StaticSurfaceInjector`, but scales the sampling
   surface with energy. This allows for much higher effective livetimes
   for low energy events in the center of the detector where veto
   techniques can be most effective.

.. cpp:class:: Floodlight

   A demonstration generator that simply illuminates its sampling surface
   from all directions with single muons.

To generate muon bundles, you must first configure a :cpp:class:`Generator` and
pass it to the provided :py:func:`GenerateBundles` segment. For example, to
generate single muons from the flux of `Hoerandel`_::

    from icecube.icetray import I3Units
    from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw
    from icecube.MuonGun.segments import GenerateBundles

    # Use Hoerandel as a template for generating muons
    model = load_model('Hoerandel5_atmod12_SIBYLL')
    # Generate only single muons, no bundles
    model.flux.max_multiplicity = 1
    # Center the sampling surface on the barycenter of IC79 strings
    surface = Cylinder(1600*I3Units.m, 800*I3Units.m, dataclasses.I3Position(31.25, 19.64, 0))
    # Draw energies from an E^-2 power law broken at 1 TeV, from 10 TeV to 10 PeV
    spectrum = OffsetPowerLaw(2, 1*I3Units.TeV, 10*I3Units.TeV, 10*I3Units.PeV)
    # Set up the generator. This gets stored in a special frame for later reference
    generator = StaticSurfaceInjector(surface, model.flux, spectrum, model.radius)

    tray.AddSegment(GenerateBundles, 'MuonGenerator', Generator=generator, NEvents=10000, GCDFile=gcd)

.. _MUPAGE: http://dx.doi.org/10.1016/j.cpc.2008.07.014
.. _Hoerandel: http://dx.doi.org/10.1016/S0927-6505(02)00198-6
