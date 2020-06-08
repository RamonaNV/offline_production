.. _ppc:

PPC (Photon Propagation Code)
=============================


This page was written by Martin Wolf with minor updates from Dima. Some of the content
is outdated, i.g., the cmake no longer makes the ppc-cli. It is still possible to compile
the stand-alone executable with "make gpu" or "make cpu" in the private/ppc/gpu directory.
Please refer to the web pages linked below for more explanation on using the stand-alone
ppc application. Much of the discussion there also applies to the icetray module.

To start using the icetray module right away, take a look at the run.py script that can
be configured to process nugen/corsika data and simulate flashers/standard candles. Use
the script "make.sh" to compile all parts of the icetray module at once.

Introduction
------------

PPC was written by Dima (Dmitry Chirkin). It computes the propagation of (Cherenkov)
photons inside the ice, in a "real time" fashion, i.e. without the use of binning and
tables (as is done in `Photonics <http://wiki.icecube.wisc.edu/index.php/Photonics>`_).
It uses the six-parameter ice model from Kurt Woschnagg [WOS2006]_. For an
explanation of the used formulas in ice.cxx see the paper by Dima [CHI2009]_.
For the Cerenkov light production a parametrization by C.H. Wiebusch [WIE1995]_
(pages 92-96) is used. 

The project's website by Dima can be found at `http://icecube.wisc.edu/~dima/work/WISC/ppc/ <http://icecube.wisc.edu/~dima/work/WISC/ppc/>`_.

There is also a `wiki page <http://wiki.icecube.wisc.edu/index.php/PPC>`_ about ppc.

PPC can be used as standalone application (ppc-cli) and as IceTray module (i3ppc).

Using PPC as standalone application
-----------------------------------

PPC as standalone application (ppc-cli) can be used for three different things:

  1. printing out the used ice model parameters
     By invoking "WFLA=wavelength ppc-cli -", e.g. "WFLA=380 ppc-cli -", with the
     wave length of the light in nanometers, ppc prints out a table
     with depth, absorption coefficient, scattering coefficient, respectively.
  2. Read in charged tracks from STDIN and produce photons according to these tracks,
     propagate them and print out hits to STDOUT when photons hit a DOM.
     By invoking "ppc-cli seed", ppc gets the seed value for the random generator
     from the command line argument. Afterwards, it reads in the track data from
     STDIN. Each line represents one track. The track format is defined by
     the F2K format [F2K2000]_ and is as follows::

       TR int int name x y z theta phi length energy time

     TR stands for track and is a character constant. The first two integer values are
     not used by ppc. The **name** column specifies the track type. Possible
     values are: "amu+", "amu-" and "amu" for muons, "delta", "brems", "epair",
     "e+", "e-" and "e" for electromagnetic cascades and "munu" and "hadr" for
     hadronic cascades. **x**, **y** and **z** are the vector components of the track's
     initial position in meters. **theta** and **phi** is the track's theta and phi angle
     in degree, respectively. **length** is the length of the track in meter.
     It is only required for muons because cascades are treated as point-like sources.
     **energy** is the track's initial energy in GeV. **time** is the track's
     initial time in nanoseconds.

     When a produced and propagated photon hits a DOM ppc prints out that hit
     in the following format::

       HIT string_no om_no time n_z

     HIT is a character constant. **string_no** is the number of the string that was hit.
     **om_no** is the number of the DOM that was hit. **time** is the time in
     nanoseconds when the DOM was hit. **n_z** is the z component of the photon's
     unit direction vector when it was hitting the DOM.
  3. Producing a flasher run by flashing one DOM.
     By invoking

       WFLA=405 ppc-cli string_no dom_no num seed

     with four arguments it produces a flasher run by flashing one specific DOM.
     **string_no** and **dom_no** is the string and DOM number, respectively.
     **num** is the number of photons that should be produced. **seed** is the
     seed value for the random number generator.

     ppc sets that DOM as source for all produced photons. It initiates the ice
     with a wavelength of 405 nm. The average theta angle of all emitted photons
     will be 90 degrees. The phi angles of the emitted photons are randomly drawn
     (0 < phi < 2PI).

     When a photon hits a DOM the same output than described in number 2 is produced.

References
----------

.. [WOS2006] Optical properties of deep glacial ice at the South Pole,
             Kurt Woschnagg et al.,
             Journal Of Geophysical Research, VOL. 111, D13203, 2006

.. [CHI2009] Study of ice transparency with IceCube flashers,
             Dmitry Chirkin,
             IceCube internal report 200911002

.. [WIE1995] Ph.D. thesis, The Detection Of Faint Light in Deep
             Underwater Neutrino Telescopes,
             C.W. Wiebusch, 1995

.. [F2K2000] F2000: A standard AMANDA offline event format,
             http://www.ifh.de/~steffenp/f2000/
