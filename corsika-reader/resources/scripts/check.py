#!/usr/bin/env python

description="""Check CORSIKA files for integrity, and print a listing of corrupt
files to stdout. The input files may be compressed using any scheme that
I3::dataio::open() can handle (e.g. gzip/bzip2/xz). A list of
"""

from optparse import OptionParser
parser = OptionParser(usage="Usage: %prog [OPTIONS] [FILES]", description=description)
parser.add_option("-v", "--verbose", default=False, action="store_true",
    help="Print names of files being checked to stderr")
parser.add_option("--detailed", default=False, action="store_true",
    help="Perform detailed check of particle blocks (much slower)")
opts, args = parser.parse_args()
if len(args) == 0:
	parser.error("You must specify at least one CORSIKA file to read!")


from icecube import icetray, dataclasses
from I3Tray import I3Tray
import sys
icetray.load('corsika-reader', False)

if not opts.verbose:
	icetray.logging.I3Logger.global_logger = icetray.I3NullLogger()
else:
	icetray.logging.set_level_for_unit('I3Tray', 'WARN')
	icetray.logging.set_level_for_unit('Python', 'INFO')

infiles = []
from glob import glob
for arg in args:
	if '*' in arg:
		infiles += sorted(glob(arg))
	else:
		infiles.append(arg)
icetray.logging.log_info('%d arguments expanded to %d files' % (len(args), len(infiles)))

def good(fname):
	try:
		tray = I3Tray()
		tray.Add('I3GSLRandomServiceFactory', Seed=1337)
		tray.Add('I3CORSIKAReader', 'reader', filenamelist=[fname], CheckIntegrity=not opts.detailed)
		tray.Add('CORSIKAResampler', 'resample',  
			OverSampling = 10,
			CylinderHeight = 1200,
			CylinderRadius = 600)
		tray.Execute()
	except RuntimeError:
		return False
	return True

baddies = 0
for f in infiles:
	if not good(f):
		print(f)
		if opts.verbose:
			sys.stderr.write('%s BAD\n' % f)
		baddies += 1
	else:
		if opts.verbose:
			sys.stderr.write('%s OKAY\n' % f)


icetray.logging.log_info('%d of %d files corrupt (%f%%)' % (baddies, len(infiles), 100*float(baddies)/len(infiles)))

if baddies > 0:
	sys.exit(1)

