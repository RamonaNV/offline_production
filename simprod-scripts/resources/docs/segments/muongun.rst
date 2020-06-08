.. _GenerateCosmicRayMuons:

GenerateCosmicRayMuons
======================

Corresponding **IceProd module**: :ref:`MuonGunGenerator`

MuonGun simulations for IceCube
-------------------------------

The ``GenerateCosmicRayMuons`` tray segment generates single atmospheric
muons (no muon bundles) with `MuonGun`_. Instead of simulating the
entire air shower, only single muon events are injected from
:math:`\theta = 0 ... 180^{\circ}` based on a cylindrical injection
surface with radius :math:`r = 800\;\rm m` and length :math:`l =
1600\;\rm m`, aligned to the z-axis, and centered at :math:`(0,0,0)\;\rm
m` in detector coordinates per default.  The in-ice muon spectrum is
approximated with a power law,

.. math::
    dP/dE_{\mu} \approx (E_{\mu} + b)^{-\gamma},

with an offset energy :math:`b = 700\;\rm GeV` and a spectral index
:math:`\gamma = 2` per default for standard IceCube simulations. The
default muon energy range is set to :math:`10\;\rm TeV ... 10\;\rm PeV`.
For each generated event, a weight is calculated according to the
*Hoerandel* model.

.. seealso::
    :ref:`MuonGun flux parametrizations` for list of available
    cosmic-ray flux parametrizations.

.. _MuonGun: ../../MuonGun/index.html


MuonGun simulations for DeepCore
--------------------------------

Low-energy simulations may set ``use_inner_cylinder=True`` in order to
switch from a static injection to an energy-dependent one: a second
smaller cylinder around DeepCore, :math:`\vec{p}_{0} =
(46.3,-34.9,-300)\;\rm m`, with radius :math:`r = 150\;\rm m` and length
:math:`l = 500\;\rm m` is created, and the injection surface is scaled
between the outer and inner cylinder w.r.t. energy.

MuonGun simulations for effective area calculations
---------------------------------------------------

The ``GenerateSingleMuons`` tray segment generates single muons with
`MuonGun`_'s :class:`icecube.MuonGun.Floodlight` injector, which
illuminates its sampling surface isotropically, and calculates the
effective area (one divided by fluence) for every generated event.


API
---

.. autofunction:: icecube.simprod.segments.GenerateCosmicRayMuons

.. autofunction:: icecube.simprod.segments.GenerateSingleMuons
