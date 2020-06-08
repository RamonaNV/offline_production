#!/usr/bin/env python

"""
Add all (potentially gigantic) histograms in a group of files.
"""

import dashi
import tables
import os, sys, operator, shutil

from optparse import OptionParser
parser = OptionParser(usage="%prog [OPTIONS] infiles outfile", description=__doc__)
parser.add_option("--blocksize", dest="blocksize", type=int, default=2048)

opts, args = parser.parse_args()
if len(args) < 2:
	parser.error("You must specify at least one output and one input file")
infiles, outfile = args[:-1], args[-1]

if os.path.exists(outfile):
	parser.error("%s already exists!" % outfile)

shutil.copy(infiles[0], outfile)

from collections import defaultdict
paths = defaultdict(list)
for fname in infiles[1:]:
	with tables.openFile(fname) as hdf:
		for group in hdf.walkNodes(where='/', classname='Group'):
			if 'ndim' in group._v_attrs: # a dashi histogram
				path = group._v_pathname
				paths[path].append(fname)

def histadd(sourceGroup, destGroup, blocksize=1):
	"""
	Add dashi histograms stored in HDF5 groups
	
	:param blocksize: operate on blocksize I/O chunks at a time
	"""
	for arr in '_h_bincontent', '_h_squaredweights':
		source = sourceGroup._v_children[arr]
		dest = destGroup._v_children[arr]
		chunksize = blocksize*reduce(operator.mul, dest.chunkshape)
		size = reduce(operator.mul, dest.shape)
		for i in range(0, size, chunksize):
			dest[i:i+chunksize] += source[i:i+chunksize]
	for prop in 'nentries', 'nans', 'nans_wgt', 'nans_sqwgt':
		destGroup._v_attrs[prop] += sourceGroup._v_attrs[prop]

with tables.openFile(outfile, 'a') as ohdf:
	for path, fnames in paths.iteritems():
		print(path)
		destGroup = ohdf.getNode(path)
		for fname in fnames:
			with tables.openFile(fname) as hdf:
				histadd(hdf.getNode(path), destGroup, opts.blocksize)
