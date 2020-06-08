from __future__ import print_function
import numpy

from icecube.photospline import splinetable, splinefitstable

from glob import glob
import re, os, sys
import copy

try:
	input = raw_input
except NameError:
	pass

from optparse import OptionParser

def stack_tables(tablist, order = 2):
	# We expect an array of (splinetable, coordinate) tuples

	bigtab = None

	for table in tablist:
		slice = table[0]
		position = table[1]

		slice.coefficients = slice.coefficients.reshape( \
		    slice.coefficients.shape + (1,))
		ndim = slice.coefficients.ndim

		if bigtab is None:
			bigtab = slice
			bigtab.knots.append([position])
			bigtab.periods.append(0)
			bigtab.order.append(order)
		else:
			bigtab.knots[ndim - 1].append(position)
			bigtab.coefficients = numpy.concatenate(
			    (bigtab.coefficients, slice.coefficients),
			    ndim - 1)

	# Shift the knots (see bsplineinterp.py)
	baseknots = bigtab.knots[ndim - 1]
	baseknots = baseknots + (numpy.max(baseknots)-numpy.min(baseknots))/(2.0*len(baseknots))*(order-1)
	interpknots = []
	for i in range (order,0,-1):
		interpknots.append(baseknots[0] - i*(baseknots[1] - baseknots[0]))
	interpknots.extend(baseknots)
	interpknots.append(interpknots[len(interpknots)-1] + (interpknots[len(interpknots)-1] - interpknots[len(interpknots)-2]))
	bigtab.knots[ndim - 1] = numpy.asarray(interpknots)

	return bigtab

def unique(seq, idfun=None):  
	"""Order-preserving uniquification, cribbed from http://www.peterbe.com/plog/uniqifiers-benchmark"""
	if idfun is None: 
		def idfun(x): return x 
	seen = {} 
	result = []
	for item in seq: 
		marker = idfun(item) 
		# in old Python versions: 
		# if seen.has_key(marker) 
		# but in new ones: 
		if marker in seen: continue 
		seen[marker] = 1 
		result.append(item)
	return result

if __name__ == "__main__":
	# Parse command line options
	parser = OptionParser(usage="%prog [options] [source dir] [fits output]")
	parser.add_option("-z", "--zstep",
	                  type="int",
	                  dest="zstep",
                      metavar="STEP",
	                  default=0,
	                  help="Increment between source depths (in meters)")
	parser.add_option("-a", "--astep",
	                  type="int",
	                  dest="astep",
                      metavar="STEP",
	                  default=0,
	                  help="Increment between source azimuth angles (in degrees)")
	parser.add_option("-f", "--filter",
	                  dest="file_type",
                      metavar="EXT",
	                  default="fits",
	                  help="File extension filter to use (e.g. '.diff.fits')")

	if len(sys.argv) < 2:
		print(sys.argv)
		sys.argv.append("-h")

	(options, args) = parser.parse_args()

	if len(args) < 1:
		print("Please supply a source directory name")
		sys.exit(0)
	if len(args) < 2:
		print("Please supply an output file name")
		sys.exit(0)

	if os.path.exists(args[1]):
		if input("File %s exists. Overwrite (y/n)? " % args[1]) == 'y':
			os.unlink(args[1])
		else:
			sys.exit()

	sourcedir = args[0] + "/"
	tables = [(i, re.match(".*_z(-?\d+)_a(\d+).*", i).groups()) for i in glob("%s*%s" % (sourcedir, options.file_type))]

	# Convert tables to a useful form and sort it

	tables = sorted([(i[0], (int(i[1][0]), int(i[1][1]))) for i in tables], key=lambda tab: tab[1])

	# Read in all the actual tables

	print('Table list acquired, reading in tables...', end=' ')
	tables = [(splinefitstable.read(i[0]), i[1]) for i in tables]
	print('done')

	depths = unique([tab[1][0] for tab in tables])
	angles = unique([tab[1][1] for tab in tables])

	extents = []
	if len(angles) > 1:
		extents.append((angles[0], angles[-1]))
	if len(depths) > 1:
		extents.append((depths[0], depths[-1]))
		

	# XXX HACK: provide full support above and below by cloning the end tables
	print("HACK: cloning tables at %.2f and %.2f" % (depths[0], depths[-1]))
	gap = depths[0] - depths[1]
	bottom = [(copy.deepcopy(tab[0]), (tab[1][0] + gap, tab[1][1])) for tab in tables if tab[1][0] == tables[0][1][0]]
	gap = depths[-1] - depths[-2]
	top = [(copy.deepcopy(tab[0]), (tab[1][0] + gap, tab[1][1])) for tab in tables if tab[1][0] == tables[-1][1][0]]
	tables = bottom + tables + top

	zpos = numpy.unique([i[1][0] for i in tables])
	if options.zstep:
		zpos = [i
		        for i in range(zpos[0],
		                       zpos[len(zpos)-1] + options.zstep,
		                       options.zstep)
		        if i in zpos
		       ]

	if len(zpos) < (zpos[len(zpos)-1] - zpos[0]) / options.zstep + 1:
		print("Error: Some depth steps are missing in table directory.")
		sys.exit(1)

	intermedtables = []

	for z in zpos:
		print('Stacking tables at z =', z)
		# Select all the tables at this z
		sublist = [tab for tab in tables if tab[1][0] == z]
		# Reformat to just one coordinate for stacking
		sublist = [(tab[0], tab[1][1]) for tab in sublist]
		if options.astep:
			sublist = [i
			           for i in sublist
			           if i[1] in range(sublist[0][1],
			                            sublist[len(sublist)-1][1] + options.astep,
			                            options.astep)
			          ]
			if len(sublist) < (sublist[len(sublist)-1][1] - sublist[0][1]) / options.astep + 1:
				print("Error: Some azimuth steps are missing in table directory.")
				print("Just stopping at z=%d" % z)
				break
				sys.exit(1)

		# extend angular range by mirroring next-to-last
		# angle bins (e.g. 170 and 10 deg) to the outside
		# (e.g. 190 and -10) so that 0 and 180 will have
		# support
		print('\tExtending angular range...', end=' ')
		lowmirror  = [(copy.deepcopy(sublist[ 1][0]), -sublist[ 1][1])]
		highmirror = [(copy.deepcopy(sublist[-2][0]),  sublist[-1][1] + (sublist[-1][1] - sublist[-2][1]))]
		sublist = lowmirror + sublist + highmirror
		print('done')

		print('\tStacking...', end=' ')
		intermedtables.append((stack_tables(sublist), z))
		print('done')

	# We no longer need to original tables
	del tables

	if len(zpos) > 1:
		finaltab = stack_tables(intermedtables)
	else:
		finaltab = intermedtables[0][0]

	try:
		targetfile = args[1]
	except IndexError:
		targetfile = os.path.normpath(os.path.join(os.path.abspath(sourcedir), '..', os.path.basename(os.path.abspath(sourcedir)) + '.fits'))

	finaltab.extents += extents

	try:
		splinefitstable.write(finaltab, targetfile)
		print("Output written to", targetfile)
	except Exception as inst:
		print(inst)


