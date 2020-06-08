.. _MuonGunGenerator:

MuonGunGenerator
----------------

The IceProd module ``MuonGunGenerator`` generates single atmospheric
muons with ``MuonGun`` and propagates them through the ice with
``PROPOSAL`` by calling the :ref:`GenerateCosmicRayMuons` and
:ref:`PropagateMuons` tray segments in its execution method,
respectively. Finally, the module writes an :cpp:class:`I3File`
containing only Q-frames; each Q-frame contains the following frame
objects:

* ``I3MCTree``
* ``MMCTrackList``
* ``RNGState``
* ``muongun_weights``


Example
-------

A basic example script; the module's parameters can be changed via the
command line::

    #!/usr/bin/env python
    # -*- coding: utf-8 -*-

    import logging

    import icecube.simprod.modules

    logging.basicConfig(level=logging.INFO)


    if __name__ == "__main__":
        mu_gun = icecube.simprod.modules.MuonGunGenerator()

        # Change some default values.
        parser = mu_gun._opt_parser
        parser.defaults["nevents"] = 100
        parser.defaults["outputfile"] = "muongun_fullchain.i3.gz"
        parser.defaults["summaryfile"] = "muongun_fullchain.xml"

        stats = {}
        mu_gun.ExecuteOpts(stats)

        for key, val in stats.iteritems():
            print "%s: %s" % (key, val)


API
---

.. autoclass:: icecube.simprod.modules.MuonGunGenerator
    :members:
    :show-inheritance:
