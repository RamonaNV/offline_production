..
.. Copyright (c) 2013
.. Claudio Kopper <claudio.kopper@icecube.wisc.edu>
.. and the IceCube Collaboration <http://www.icecube.wisc.edu>
..
.. Permission to use, copy, modify, and/or distribute this software for any
.. purpose with or without fee is hereby granted, provided that the above
.. copyright notice and this permission notice appear in all copies.
..
.. THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
.. WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
.. MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
.. SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
.. WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
.. OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
.. CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
..
..
.. $Id: index.rst 110958 2013-09-19 20:55:39Z claudio.kopper $
..
.. @file index.rst
.. @version $Revision: 110958 $
.. @date $Date: 2013-09-19 15:55:39 -0500 (Thu, 19 Sep 2013) $
.. @author Claudio Kopper
..

.. highlight:: python

.. _fullSimulation-main:

============================
How to use fullSimulation.py
============================

This is mostly copy&pasted from an email written to Donglian a while ago
and could probably use an overhaul.

Prerequisites
~~~~~~~~~~~~~

In order to make your own nu_tau (or nu_anything) sample,
you just need the most recent release candidate of "icesim4".
You can find that here: http://code.icecube.wisc.edu/svn/meta-projects/simulation/candidates/TODO

We try to keep simulation trunk stable, so you can also try that if you run into problems:
http://code.icecube.wisc.edu/svn/meta-projects/simulation/trunk

The scripts use clsim, so you will need a system with OpenCL installed
(it does not need a GPU, so any Madison system should work just fine.
OpenCL on Macs also just works.).
It should tell you during the "cmake" stage of clsim if OpenCL can or cannot be found.

As always, just check it out, compile it and you should be old to go.
The simulation scripts are included; they are in "simprod-scripts" as tray segments.
If you want to have a look at the individual parts, they are in ``$I3_SRC/simprod-scripts/python/segments/`` .

The main script putting them all together is here: ``$I3_SRC/simprod-scripts/resources/examples/fullSimulation.py``

You should be able to run that directly to generate your sample
(you need to change the first couple of lines to point it to some
temporary scratch space and to the spline tables if you want to do hybrid simulation
(Hybrid simulation cannot handle tilted ice (all layers will be flat) and it cannot handle ice anisotropy at all,
so SpiceLea will not work. In case you choose not to use it, you can just ignore the cascade table path setting).
The current paths are set up for Madison systems.


Some Example Options
~~~~~~~~~~~~~~~~~~~~

This is an example command line used to generate a couple of nu_tau events::

  ./fullSimulation.py -n 10 --seed=43682 --datasetnumber=1 --runnumber=1 --no-hybrid --icemodel=SpiceLea --detector=IC86 --unshadowed-fraction=0.99 --flavor=NuTau --outfile=taus.i3 --skip-calibration --from-energy=1000 --to-energy=10000000 --include-gcd-in-outfile

I'll go through what these options mean: (you can use ``fullSimulation.py --help`` to get some information, too)

* ``-n 10``
  Simulate 10 events (the output file might contain fewer events because not everything will trigger).

* ``--seed=43682``
  A random number generator seed

* ``--datasetnumber=1 --runnumber=1``
  A combination of "dataset" and "run" numbers. You can keep the same seed for all of them, each one will give a distinct independent set of events.
  (So you can think of the combination of seed, dataset and run as the actual random number seed. Or, put another way, there is no need to have different 
  seeds for different runs.)

* ``--no-hybrid``
  Do not use hybrid simulation mode, i.e. propagate everything using direct simulation. There will be no photon tables and things like ice anisotropy/SpiceLea will work.

* ``--icemodel=SpiceLea``
  The other two models are "Spice1" and "SpiceMie". Including things like WHAM would be trivial if you need it.

* ``--detector=IC86``
  This will select a GCD file from ``$I3_TESTDATA`` automatically. Currently works for IC86 and IC79.

* ``--unshadowed-fraction=0.99``
  This is the "DOMEfficiency", currently named like this for compatibility with other tools that use the same name.

* ``--flavor=NuTau``
  You can set this also to "NuE" and "NuMu".

* ``--outfile=...``
  The name of your final .i3 file you want to generate.

* ``--from-energy=1000 --to-energy=10000000``
  The energy range in GeV. Currently not implemented for MuonGun which uses fixed energy ranges. Should be fixed soon.

* ``--include-gcd-in-outfile``
  Use this option if you want to generate output files with GCD frames in them. It makes them much easier to use, but of course wastes some space..

-----------

And here are some more options you might want to use:

* ``--outfile-gcd=gcd_file.i3``
  You can use this instead of "--include-gcd-in-outfile" to write GCD frames to a separate file instead of the output file. You can, for example, include this in "run 0" only to make one GCD file and then not use it for all the other runs if you make more than one.

* ``--use-gpu``
  If you want to run this on GPUs (you can do it interactively on "cobaltgpu" and "cobaltamdgpu"). In that case you should also set "--max-parallel-events=100" or some number between 10 and 100. This will speed things up on GPUs considerably but use more system memory. I think "10" should be okay for high-energy simulations. For low energies you would use something more like 100. and for PINGU it would be 1000. It basically bundles events together on the GPU to use less bandwidth..


Summary
~~~~~~~

This should be a tool to generate your own neutrino MC files with whatever settings you need.
You can send jobs to the cluster (although simulation will be very slow on CPUs)
or just make a couple of files on GPU machines.
I would try cobaltgpu and/or cobaltamdgpu for testing - if it works, we have a lot of GPU nodes on npx4 now.
The cobaltgpu machines even have their own local condor queues, so you can submit jobs
to the GPUs and only four of them will ever run in parallel (they both have four GPUs).
There's an example submit script for cobaltgpu here: https://wiki.icecube.wisc.edu/index.php/CobaltGPU
You get an "interactive" job when you follow these instructions,
i.e. it just opens a new shell with your environment set to make it only use a single GPU.
(You might not realize you just started a job.. ;-) Typing ``exit`` terminates it again.)

