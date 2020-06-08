..
.. copyright  (C) 2012
.. The Icecube Collaboration
..
.. $Id: index.rst 141130 2016-01-21 20:40:46Z david.schultz $
..
.. @version $Revision: 141130 $
.. @date $LastChangedDate: 2016-01-21 14:40:46 -0600 (Thu, 21 Jan 2016) $
.. @author Nathan Whitehorn <nwhitehorn@physics.wisc.edu> $LastChangedBy: david.schultz $

.. _corsika-resampler:

CORSIKAResampler
================

CorsikaResampler Module
-----------------------

CorsikaResampler provides a module to resample events read by I3CORSIKAReader. It then
Randomly shifts the shower core around the surface of a cylinder with configurable dimensions.
Each simulated shower is resampled N times and emitted as a single
DAQ frame containing an I3MCTree with the shower particles, as well as an
object containing weight information and a header with the shower number.

Arguments
^^^^^^^^^

1. OverSampling (default: ``1``)

2. CylinderHeight (default: ``1600`` m)

	Height of vertical cylinder centered at (0,0,0) in detector coordinates
	through which all simulated shower cores are made to pass.

3. CylinderRadius (default: ``800`` m)

	Radius of vertical cylinder centered at (0,0,0) in detector coordinates
	through which all simulated shower cores are made to pass.

6. ZenithBias (default: ``True``)

	Set to true if CORSIKA was compiled with the VOLUMECORR option and false
	if the VOLUMEDET option was used. The default zenith bias (proportional
	to sin(theta) for surface detectors) is not yet supported.

