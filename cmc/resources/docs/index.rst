..
.. Copyright (c) 2020
.. Bernhard Voigt <bernhard.voigt@desy.de>
.. Sebastian Panknin <panknin@physik.hu-berlin.de>
.. Alex Olivas <aolivas@umd.edu>
.. Juan Carlos Diaz-Velez <juancarlos.diazvelez@icecube.wisc.edu>
.. Justin Lafranchi <jll1062@psu.edu>
.. Brian Clark <brianclark@icecube.wisc.edu>
..
.. Permission to use, copy, modify, and/ordistribute this software for any
.. purpose with or without fee is hereby granted, provided that the above
.. copyright notice and this permission notice appear in all copies.
..
.. THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
.. WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
.. MERCHANTABILIITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
.. SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
.. WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
.. OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
.. CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
..
..
.. @file index.rst
.. @version $LastChangedRevision$
.. @date $Date$
.. @author Brian Clark

.. _cmc:

cmc - Cascade Monte Carlo
=========================

The cmc (Cascade Monte Carlo) module simulates the development of
electro-magnetic showers in the ice.

This module should be used after the neutrino-generator in case of
electron neutrion samples. If you want to treat cascades from muon
interactions, put into the MC tree by mmc, you need to run this module
after PROPOSAL.  But if you want to look at hadronic cascades and use
this module for generating muons, you have to run PROPOSAL after it.

The original documentation for his package was wiki based here: http://wiki.icecube.wisc.edu/index.php/Cascade_Simulation. That documentation is no longer maintained, and is kept only for historical purposes.

.. toctree::
   :maxdepth: 3
   
   release_notes
   cmc_code
   physics_overview
