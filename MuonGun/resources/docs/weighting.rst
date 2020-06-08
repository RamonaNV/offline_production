Weighting
=========

In order to weight generated muon bundles to a flux, you need to know both the
distribution of events in the flux model (parameterized in the tables provided
with MuonGun) and the distribution of events that you generated (calculated by
the Generator, stored in a special "S" frame at the beginning of every
generated file). Given those, you can calculate weights either within IceTray
using the :cpp:class:`WeightCalculatorModule` I3Module or from a standalone
Python script using the :cpp:class:`WeightCalculator` class.

First, you should collect the generators for all of the files you plan to use::

    def harvest_generators(infiles):
        """
        Harvest serialized generator configurations from a set of I3 files.
        """
        from icecube.icetray.i3logging import log_info as log
        generator = None
        for fname in infiles:
            f = dataio.I3File(fname)
            fr = f.pop_frame(icetray.I3Frame.Stream('S'))
            f.close()
            if fr is not None:
                for k in fr.keys():
                    v = fr[k]
                    if isinstance(v, MuonGun.GenerationProbability):
                        log('%s: found "%s" (%s)' % (fname, k, type(v).__name__), unit="MuonGun")
                        if generator is None:
                            generator = v
                        else:
                            generator += v
        return generator

The generators can simply be added together, or multiplied by an integer to
represent a larger number of identically-configured generators.

Weighting to a cosmic ray flux model
------------------------------------

You can pass
this combined generator to :cpp:class:`WeightCalculatorModule` to calculate a
weight appropriate for the combined set of files::

    model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
    generator = harvest_generators(infiles)
    tray.AddModule('I3MuonGun::WeightCalculatorModule', 'MuonWeight', Model=model,
        Generator=generator)

This will put an I3Double called "MuonWeight" into the frame that represents a
weight in events per second. Alternatively, you can use the provided
:ref:`tableio-main` converter to write the parameters needed for the weight
calculation to a table::

    from icecube.hdfwriter import I3HDFWriter
    tray.AddSegment(I3HDFWriter, 'scribe',
        Output=outfile,
        Keys=[dict(key='I3MCTree', name='BundleParameters',
                 converter=MuonGun.converters.MuonBundleConverter(1, generator.surface))],
        Types=[],
        SubEventStreams=['nullsplit'],
    )

and then use the standalone :cpp:class:`WeightCalculator` class to calculate a
weight::

    model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
    generator = harvest_generators(infiles)
    weighter = MuonGun.WeightCalculator(generator.surface, model, generator)

    with tables.openFile(outfile) as hdf:
        axis = hdf.root.MCPrimary.read()
        bundle = hdf.root.BundleParameters.read()
        weights = weighter(axis['x'], axis['y'], axis['z'], axis['zenith'], axis['azimuth'],
            bundle['multiplicity'], bundle['energy'], bundle['radius'])

.. note:: The weighter will only be able to accept Numpy arrays if you have `boost::numpy`_ installed. If you do not have `boost::numpy`_ it will simply be exposed as a scalar function.

.. _`boost::numpy`: https://github.com/martwo/BoostNumpy/

Calculating a muon effective area
---------------------------------

To calculate a muon effective area for single muons, simply sum up the inverse
of the generated fluences for each event, e.g.::

	mctree = frame['I3MCTree']
	primary = mctree.primaries[0]
	muon = mctree.get_daughters(primary)[0]
	bundle = BundleConfiguration([BundleEntry(0, muon.energy)])
	area_weight = 1./generator.generated_events(primary, bundle)

``area_weight`` has units of :math:`GeV \, m^{2} \, sr`, and is analogous to
NeutrinoGenerator's ``OneWeight``. To obtain the effective area in units of
:math:`m^{2}` as a function of muon energy and direction, fill ``area_weight``
into a histogram and divide each bin by its width in energy and solid angle.
