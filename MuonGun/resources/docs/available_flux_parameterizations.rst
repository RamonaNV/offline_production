.. _MuonGun flux parametrizations:

Available flux parameterizations
================================

The following models are distributed with MuonGun, and are valid arguments for
:py:func:`icecube.MuonGun.load_model`:

.. table:: Primary cosmic-ray flux parameterizations

   +--------------------------------+---------+-----------+----------------+
   | Model string                   |CR flux  |Atmosphere | Hadronic model |
   +================================+=========+===========+================+
   | Hoerandel5_atmod12_SIBYLL      |Hoerandel|12 (winter)|SIBYLL          |
   +--------------------------------+---------+           |                |
   | GaisserH4a_atmod12_SIBYLL      |Gaisser  |           |                |
   +--------------------------------+         |           +----------------+
   | GaisserH4a_atmod12_DPMJET      |         |           |DPMJET (conv.)  |
   +--------------------------------+         |           +----------------+
   | GaisserH4a_atmod12_DPMJET-C    |         |           |DPMJET (prompt) |
   +--------------------------------+---------+-----------+----------------+

In addition to the cosmic-ray flux models there are also two 'pseudofluxes' that
parameterize the output of dCORSIKA in two configurations used to generate most
of the IC79 penetrating muon simulation.

.. table:: CORSIKA output parameterizations

   +--------------------------------------+
   | Model string                         |
   +======================================+
   | Standard5Comp_atmod12_SIBYLL         |
   +--------------------------------------+
   | CascadeOptimized5Comp_atmod12_SIBYLL |
   +--------------------------------------+
