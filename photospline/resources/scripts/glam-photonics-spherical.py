#!/usr/bin/env python

from icecube.photospline import splinefitstable
from optparse import OptionParser

from icecube.photospline.photonics import *

import sys
import os
import numpy

try:
	input = raw_input
except NameError:
	pass

# Hard-coded params

#nknots =[17, 6, 12, 25] # [r, phi, z, t]  For Nathan/Jakob's binning

# Parse arguments

usage = "usage: %prog [options] table.pt [output.fits]"
optparser = OptionParser(usage=usage)
optparser.add_option("-r", "--rknots", dest="rknots", type="int",
             help="number of knots in radial dimension")
optparser.add_option("-f", "--fknots", dest="fknots", type="int",
             help="number of knots in angular dimension")
optparser.add_option("-z", "--zknots", dest="zknots", type="int",
             help="number of knots in cos(zenith)")
optparser.add_option("-t", "--tknots", dest="tknots", type="int",
             help="number of knots in time dimension")
optparser.add_option("-s", "--smooth", dest="smooth", type="float",
             help="smoothness coefficient", default=1e-6)
optparser.add_option("--prob", dest="prob", action="store_true",
             help="Fit only the normalized CDFs", default=False)
optparser.add_option("--abs", dest="abs", action="store_true",
             help="Fit only the total amplitude in each cell", default=False)
optparser.add_option("--force", dest="force", action="store_true",
             help="Overwrite existing fits files", default=False)
optparser.add_option("--ice-bottom", dest="ice_bottom", type="float",
             help="Lower boundary of ice properties. Any table cells below this\
             depth will be weighted with zero, as they contain no data.", default=-820)
optparser.add_option("--ice-top", dest="ice_top", type="float",
             help="Upper boundary of ice properties.", default=820)
(opts, args) = optparser.parse_args()

if len(args) < 1:
	optparser.print_usage()
	sys.exit(1)

if (not os.path.exists(args[0])):
	optparser.error("Input table %s doesn't exist!" % args[0])

if os.stat(args[0]).st_size == 0:
	optparser.error("Input table %s has zero size! Did photomc finish properly?" % args[0])

# by default, do both fits
if not opts.prob and not opts.abs:
	opts.prob = opts.abs = True

def check_exists(outputfile):
    if os.path.exists(outputfile):
        if opts.force or input("File %s exists. Overwrite? (y/n)" % outputfile) == 'y':
            os.unlink(outputfile)
        else:
            sys.exit()

def default_path(input):
	pth = os.path.basename(input)
	return pth+'.abs.pspl.fits',pth+'.prob.pspl.fits'

if len(args) < 2:
    abs_outputfile, prob_outputfile = default_path(args[0])
else:
	if opts.prob and opts.abs:
		# Output must be base name
		abs_outputfile, prob_outputfile = default_path(args[1])
	else:
		# Name whichever the exact name
		abs_outputfile = prob_outputfile = args[1]

smooth = opts.smooth

# Real code
from icecube.photospline import spglam as glam

# Load original table. If it fails, try another format.
try:
	table = Table(args[0], mmap=False)
	geo = table.ph_header.geo
except:
	table = FITSTable.load(args[0])
	geo = table.header['geometry']

table.normalize()

# check for a sane normalization
eff = Efficiency.RECEIVER | Efficiency.WAVELENGTH | Efficiency.N_PHOTON
int2bin = lambda num, count: "".join([str((num >> y) & 1) for y in range(count-1, -1, -1)])
if (table.header['efficiency'] != eff):
	err = "Unknown normalization %s (expected %s)" % (
	    int2bin(table.header['efficiency'], 9), int2bin(eff, 9))
	raise ValueError(err)

if (table.header['geometry'] is not Geometry.SPHERICAL):
	raise ValueError("This table does not have spherical geometry")

# Extents for a particular set of spherically-binned tables produced in 2010
extents = [(0.0, 600.0), (0.0, 180.0), (-1.0, 1.0), (0.0, 7000.)] # r, phi, cos(polar), t
ngroup = table.header['n_group']

def construct_knots(nknots = None):
	if nknots is None:
		nknots = [15, 6, 15]
		if table.ndim > 3:
			nknots.append(25) # [t]
	
	if opts.rknots:
	    nknots[0] = opts.rknots
	if opts.fknots:
	    nknots[1] = opts.fknots
	if opts.zknots:
	    nknots[2] = opts.zknots
	if opts.tknots and table.ndim > 3:
	    nknots[3] = opts.tknots
	
	print("Core knots:", nknots)
	
	radial_extent = extents[0][1]
	coreknots = [None]*4
	
	# It's tempting to use some version of the bin centers as knot
	# positions, but this should be avoided. Data points exactly at the
	# knot locations are not fully supported, leading to genuine wierdness
	# in the fit.
	coreknots[0] = numpy.linspace(extents[0][0], radial_extent**(1./2),
	    nknots[0])**2
	coreknots[1] = numpy.linspace(extents[1][0], extents[1][1], nknots[1])
	coreknots[2] = numpy.linspace(extents[2][0], extents[2][1], nknots[2])
	
	# We're fitting the CDF in time, so we need tightly-spaced knots at
	# early times to be able to represent the potentially steep slope.
	coreknots[3] = numpy.logspace(0, numpy.log10(extents[3][1]), nknots[3])
	coreknots[3] = numpy.concatenate(([0], coreknots[3]))
	
	# Now append the extra knots off both ends of the axis in order to
	# provide full support at the boundaries
	
	rknots     = numpy.append(numpy.append([-1, -0.5, -0.1], coreknots[0]),
	                          100*numpy.arange(1,3) + radial_extent)
	endgap = [coreknots[1][1]-coreknots[1][0],
	    coreknots[1][-1]-coreknots[1][-2]]
	thetaknots = numpy.concatenate((coreknots[1][0] -
	    endgap[0]*numpy.arange(2,0,-1), coreknots[1], coreknots[1][-1] +
	    endgap[1]*numpy.arange(1,3)))
	# NB: we want -1 and 1 to be fully supported.
	endgap = [coreknots[2][1]-coreknots[2][0],
	    coreknots[2][-1]-coreknots[2][-2]]
	zknots = numpy.concatenate((coreknots[2][0] - 
	    endgap[0]*numpy.arange(2,0,-1), coreknots[2], coreknots[2][-1] +
	    endgap[1]*numpy.arange(1,3)))
	
	# NB: we can get away with partial support in time, since we know that
	# F(0) is identically zero.
	tknots = numpy.concatenate((coreknots[3], coreknots[3][-1] +
	    100*numpy.arange(1,4)))
	
	print('knots:')
	print(rknots)
	print(thetaknots)
	print(zknots)
	print(tknots)

	return [rknots, thetaknots, zknots, tknots]

def spline_spec(ndim):
   if ndim > 3:
       order = [2,2,2,3]        # Quadric splines for t to get smooth derivatives
       penalties = {2:[smooth]*3 + [0], # penalize curvature in rho,z,phi
                    3:[0]*3 + [smooth]} # order 3 in time CDF => order 2 in time PDF
       knots = construct_knots([20, 8, 15, 25])
   else:
       order = [2,2,2]    # Quadric splines to get smooth derivatives
       penalties = {2:[smooth]*3}    # Penalize curvature 
       knots = construct_knots([25, 25, 25, 25])[:3]
   return order, penalties, knots

# Take cumulative sum to get the CDF, and adjust fit points to be
# the right edges of the time bins, where the CDF is measured.
table.values[:] = numpy.cumsum(table.values, axis=3)
table.bin_centers[3] += table.bin_widths[3]/2.

print("Loaded histogram with dimensions ", table.shape)

norm = table.values[:,:,:,-1]

# Rescale all axes to have a maximum value of ~ 10
def rescale_axes(knots, bin_centers, bin_widths):
	axis_scale = []
	for i in range(0,len(bin_centers)):
		scale = 2**numpy.floor(numpy.log(numpy.max(bin_centers[i])/10.) /
		    numpy.log(2))
		axis_scale.append(scale)
		bin_centers[i] /= scale
		knots[i] /= scale
		bin_widths[i] /= scale
	return axis_scale

if opts.abs:
	z = numpy.log(norm)

	# add some numerical stability sauce
	w = 1000*numpy.ones(norm.shape)
	w[numpy.logical_not(numpy.isfinite(z))] = 0
	z[numpy.logical_not(numpy.isfinite(z))] = 0
	
	# XXX HACK: don't believe anything that happens outside the
	#           tracking volume of the table
	#scalp(table, w, low=opts.ice_bottom, high=opts.ice_top)

	order, penalties, knots = spline_spec(3)
	bin_centers = [b.copy() for b in table.bin_centers[:3]]
	bin_widths = [b.copy() for b in table.bin_widths[:3]]
	axis_scale = rescale_axes(knots, bin_centers, bin_widths)

	print('Number of knots used: ',[len(a) for a in knots])
	print("Beginning spline fit for abs table...")
	spline = glam.fit(z,w,bin_centers,knots,order,smooth,penalties=penalties)
	spline.geometry = Geometry.SPHERICAL
	spline.extents = extents[:3]
	spline.ngroup = table.header['n_group']

	print("Saving table to %s..." % abs_outputfile)
	spline.knots = [spline.knots[i] * axis_scale[i] for i
			    in range(0, len(spline.knots))]
	check_exists(abs_outputfile)
	splinefitstable.write(spline, abs_outputfile)

	# clean up
	del(w,z,bin_centers,bin_widths,order,penalties,knots,spline)

if opts.prob:
	z = table.values / norm.reshape(norm.shape + (1,))
	# XXX HACK: ignore weights for normalized timing
	w = 1000*numpy.ones(table.values.shape)

	# XXX HACK: don't believe anything that happens outside the
	#           tracking volume of the table
	#scalp(table, w, low=opts.ice_bottom, high=opts.ice_top)

	order, penalties, knots = spline_spec(4)
	bin_centers = [b.copy() for b in table.bin_centers]
	bin_widths = [b.copy() for b in table.bin_widths]
	axis_scale = rescale_axes(knots, bin_centers, bin_widths)

	# go ahead and remove the table from memory
	del(table, norm)

	print('Number of knots used: ',[len(a) for a in knots])
	print("Beginning spline fit for timing table...")
	spline = glam.fit(z,w,bin_centers,knots,order,smooth,penalties=penalties,monodim=3)
	spline.geometry = Geometry.SPHERICAL
	spline.extents = extents
	spline.ngroup = ngroup

	print("Saving table to %s..." % prob_outputfile)
	spline.knots = [spline.knots[i] * axis_scale[i] for i
			    in range(0, len(spline.knots))]
	check_exists(prob_outputfile)
	splinefitstable.write(spline, prob_outputfile)

	# clean up
	del(w,z,bin_centers,bin_widths,order,penalties,knots,spline)


# smoothed = glam.grideval(spline, table.bin_centers)
# resid = (smoothed - table.values)[table.weights != 0]
# fracresid = ((smoothed - table.values)/table.values)[table.weights != 0]
# 
# 
# print "Fit Statistics:"
# print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(resid))
# print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean(resid**2))
# print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs(fracresid))
# print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs(fracresid))

