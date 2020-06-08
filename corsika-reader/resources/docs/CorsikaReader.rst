..
.. copyright  (C) 2012
.. The Icecube Collaboration
..
.. $Id: index.rst 141130 2016-01-21 20:40:46Z david.schultz $
..
.. @version $Revision: 141130 $
.. @date $LastChangedDate: 2016-01-21 14:40:46 -0600 (Thu, 21 Jan 2016) $
.. @author Nathan Whitehorn <nwhitehorn@physics.wisc.edu> $LastChangedBy: david.schultz $

.. _corsika-reader:

CORSIKAReader
=============

CorsikaReader Module
--------------------

I3CorsikaReader provides a module to read output binary files from the CORSIKA
air-shower simulation package. Each simulated shower is emitted as a single
DAQ frame containing an I3MCTree with the shower particles, as well as an
object containing weight information and a header with the shower number.

Arguments
^^^^^^^^^

1. FilenameList (default: ``[]``)

	Array of paths to CORSIKA DAT files to read. Will be emitted in order.
	Each CORSIKA run (from the CORSIKA run configuration) will be emitted
	with a separate run ID in the frames.

2. Prefix (default: ``None``)

	Path to I3 file (e.g. GCD information) to emit before the CORSIKA
	events.

3. NEvents (default: ``0``)

	Number of CORSIKA showers per input file. This is required only on
	older CORSIKA versions that do not write this information into the
	run header. On such CORSIKA files, it will be checked at run end and
	an error emitted if incorrect. This parameter is ignored if set in
	the file.
