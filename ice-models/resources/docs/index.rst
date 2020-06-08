.. _ice-models:

Ice Models
==========

.. toctree::
   :maxdepth: 1
   
   release_notes

Introduction
------------

The project contains ice models that can be used by both ppc and clsim.


Data Files
----------

The minimum set of files consists of:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cfg.txt
must have 4 lines:

* 5     # over-R: DOM radius "oversize" scaling factor
* 1.0   # overall DOM efficiency correction
* 0.35  # 0=HG; 1=SAM (scattering function shape parameter)
* 0.9   # g=<cos(theta)>

icemodel.dat
must have 2 or more lines containing 4 entries per line:

* depth of the center of the layer in meters
* effective scattering coefficient b_e(400)
* absorption from dust a_dust(400)
* temperature delta (vs. value at 1730 m) - as defined in section 4 of the `SPICE paper <http://arxiv.org/abs/1301.5361>`_

icemodel.par
4 or 6 lines, two values per line,
second value is uncertainty and is ignorred:

* alpha: absorption wavelength dependence exponent
* kappa: scatterint wavelength dependence exponent
* A: A and B define pure ice absorption
* B: as A*exp(-B/wavelength in nm)
* D: if present, the absorption coefficient is calculated as
* E: a_dust(400)=(D*[3rd element in a line]+E)*400^-kappa


Optionally more files can be added and existing files extended:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tilt.dat and tilt.par - these define ice layer tilt. The direction of the tilt is hard-coded in ppc at 225 degrees SW. These files are normally provided with particular ice models. In these models tilt can be disabled by removing these two files.

Anizotropy is enabled by adding the following lines to the file cfg.txt:

* 130   # direction of major anisotropy axis
* -0.106 # magnitude of major anisotropy coefficient k1
* 0.053  # magnitude of minor anisotropy coefficient k2

Additionally if your code supports hole ice (currently CUDA version of ppc does, if recompiled with #define HOLE enabled), the hole ice is defined by the following lines at the end of the file cfg.txt:

* 0.5   # hole ice radius in units of [DOM radius]
* 0.5   # hole ice effective scattering length [m]
* 100   # hole ice absorption length [m]
* 0.5   # hole ice 0=HG; 1=SAM
* 0.9   # hole ice g=<cos(theta)>

For use with ppc additional files must be added to the ice configuration directory:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* wv.dat: wavelength-tabulated DOM acceptance convolved with the emission spectrum (Cherenkov or flasher LED)
* as.dat: parameters of the angular sensitivity polynomial expansion. To record all arriving photons without applying the angular sensitivity curve, create the file as.dat with only one element set to 1 (with echo 1 > as.dat).
* rnd.txt: table of random number multipliers for the multiply-with-carry random number generator


