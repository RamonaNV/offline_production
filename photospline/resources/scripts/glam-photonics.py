from icecube.photospline import splinefitstable
from optparse import OptionParser

from icecube.photospline.photonics import *

try:
	input = raw_input
except NameError:
	pass

import sys
import os
import numpy

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
             help="number of knots in longitudinal dimension")
optparser.add_option("-t", "--tknots", dest="tknots", type="int",
             help="number of knots in time dimension")
optparser.add_option("-s", "--smooth", dest="smooth", type="float",
             help="smoothness coefficient", default=1e-6)
optparser.add_option("--prob", dest="prob", action="store_true",
             help="Fit only the normalized CDFs", default=False)
optparser.add_option("--abs", dest="abs", action="store_true",
             help="Fit only the total amplitude in each cell", default=False)
optparser.add_option("--ice-bottom", dest="ice_bottom", type="float",
             help="Lower boundary of ice properties. Any table cells below this\
             depth will be weighted with zero, as they contain no data.", default=-820)
optparser.add_option("--ice-top", dest="ice_top", type="float",
             help="Upper boundary of ice properties. Any table cells above this\
             depth will be weighted with zero, as they contain no data.", default=820)
(opts, args) = optparser.parse_args()
if len(args) < 1:
	print(usage)
	sys.exit(1)

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
    abs_outputfile, prob_outputfile = default_path(args[1])

if opts.prob: check_exists(prob_outputfile)
if opts.abs: check_exists(abs_outputfile)

smooth = opts.smooth

# Real code
from icecube.photospline import spglam as glam

table = Table(args[0])

table.convert_to_level1()

# Photonics stores a bitmask that gives the kinds of normalizations
# that have been applied to the table cells in the 'efficiency' field.
# NB: We want dP, not dP/dt
if (Efficiency.DIFFERENTIAL & table.header['efficiency']):
	raise ValueError("This appears to be a dP/dt table. Don't do that, okay?")
if (not Efficiency.N_PHOTON & table.header['efficiency']):
	raise ValueError("This table does not appear to be normalized.")

nknots = [15, 6, 25] # rho, phi, z
if table.ndim > 3:
    nknots.append(20) # [t]

if opts.rknots:
    nknots[0] = opts.rknots
if opts.fknots:
    nknots[1] = opts.fknots
if opts.zknots:
    nknots[2] = opts.zknots
if opts.tknots and table.ndim > 3:
    nknots[3] = opts.tknots

print("Core knots:", nknots)

radial_extent = 600
length_extent = 500

coreknots = [None]*4

# It's tempting to use some version of the bin centers as knot positions,
# but this should be avoided. Data points exactly at the knot locations are 
# not fully supported, leading to genuine wierdness in the fit.
coreknots[0] = numpy.linspace(0, numpy.sqrt(radial_extent), nknots[0])**2
coreknots[0] = numpy.concatenate(([0], numpy.logspace(-1,
    numpy.log10(radial_extent), nknots[0]-1)))
coreknots[1] = numpy.linspace(0, 180, nknots[1])
# space 1/3 of the knots quadratically behind the source, 
# where everything is diffuse, and the remainder in front
# with logarithmic spacing
backerds = int(nknots[2]/3.0)
coreknots[2] = numpy.concatenate((
    -(numpy.linspace(1, numpy.sqrt(length_extent), backerds)**2)[::-1],
    numpy.logspace(0, numpy.log10(length_extent), nknots[2]-backerds)
    ))

# We're fitting the CDF in time, so we need tightly-spaced knots at
# early times to be able to represent the potentially steep slope.
# XXX: we assume t_max == 7000 ns
coreknots[3] = numpy.logspace(-1, numpy.log10(7000), nknots[3])
coreknots[3] = numpy.concatenate(([0], coreknots[3]))

# Now append the extra knots off both ends of the axis in order to provide
# full support at the boundaries

rknots     = numpy.append(numpy.append([-1, -0.5, -0.1], coreknots[0]),
                          100*numpy.arange(1,3) + radial_extent)
endgap = [coreknots[1][1]-coreknots[1][0], coreknots[1][-1]-coreknots[1][-2]]
thetaknots = numpy.concatenate((coreknots[1][0] - endgap[0]*numpy.arange(2,0,-1),
    coreknots[1], coreknots[1][-1] + endgap[1]*numpy.arange(1,3)))
# NB: we want -1 and 1 to be fully supported.
endgap = [coreknots[2][1]-coreknots[2][0], coreknots[2][-1]-coreknots[2][-2]]
zknots = numpy.concatenate((coreknots[2][0] - endgap[0]*numpy.arange(2,0,-1),
    coreknots[2], coreknots[2][-1] + endgap[1]*numpy.arange(1,3)))

# NB: we can get away with partial support in time, since we know that
# F(0) is identically zero.
tknots = numpy.concatenate((coreknots[3], 7000 + 100*numpy.arange(1,4)))

print('knots:')
print(rknots)
print(thetaknots)
print(zknots)
print(tknots)

def spline_spec(ndim):
   if ndim > 3:
       order = [2,2,2,3]        # Quadric splines for t to get smooth derivatives
       penalties = {2:[smooth]*3 + [0], # penalize curvature in rho,z,phi
                    3:[0]*3 + [smooth]} # order 3 in time CDF => order 2 in time PDF
       knots = [rknots, thetaknots, zknots, tknots]
   else:
       order = [2,2,2]    # Quadric splines to get smooth derivatives
       penalties = {2:[smooth]*3}    # Penalize curvature 
       knots = [rknots, thetaknots, zknots]
   return order, penalties, knots

# Take cumulative sum to get the CDF, and adjust fit points to be
# the right edges of the time bins, where the CDF is measured.
table.values = numpy.cumsum(table.values, axis=3)
table.bin_centers[3] += table.bin_widths[3]/2.

print("Loaded histogram with dimensions ", table.shape)

norm = table.values[:,:,:,-1]

# Rescale all axes to have a maximum value of ~ 10
axis_scale = []
knots = [rknots, thetaknots, zknots, tknots]
for i in range(0,len(table.bin_centers)):
	scale = 2**numpy.floor(numpy.log(numpy.max(table.bin_centers[i])/10.) /
	    numpy.log(2))
	axis_scale.append(scale)
	table.bin_centers[i] /= scale
	knots[i] /= scale
	table.bin_widths[i] /= scale

if opts.abs:
	z = numpy.log(norm)

	# add some numerical stability sauce
	w = 1000*numpy.ones(norm.shape)

	# Zero (and remove from fit) table cells with non-finite values
	# (e.g. photon count was zero, and we have log(0) at this point)
	w[numpy.logical_not(numpy.isfinite(z))] = 0
	z[numpy.logical_not(numpy.isfinite(z))] = 0

	# XXX HACK: don't believe anything that happens outside the
	#           tracking volume of the table
	#scalp(table, w, low=opts.ice_bottom, high=opts.ice_top)
	# XXX HACK: don't believe anything in the first 3 radial bins
	#w[:3,:,:] = 0

	order, penalties, knots = spline_spec(3)

	print('Number of knots used: ',[len(a) for a in knots])
	print("Beginning spline fit for abs table...")
	spline = glam.fit(z,w,table.bin_centers[:3],knots,order,smooth,penalties=penalties)

	print("Saving table to %s..." % abs_outputfile)
	spline.knots = [spline.knots[i] * axis_scale[i] for i
			    in range(0, len(spline.knots))]
	splinefitstable.write(spline, abs_outputfile)

	# clean up
	del(w,z,order,penalties,knots,spline)

if opts.prob:
	z = table.values / norm.reshape(norm.shape + (1,))
	# Same sauce as above.
	w = 1000*numpy.ones(table.weights.shape)
	w[numpy.logical_not(numpy.isfinite(z))] = 0
	z[numpy.logical_not(numpy.isfinite(z))] = 0
	order, penalties, knots = spline_spec(4)

	centers = table.bin_centers

	# XXX HACK: don't believe anything that happens outside the
	#           tracking volume of the table
	#scalp(table, w, low=opts.ice_bottom, high=opts.ice_top)
	# XXX HACK: also, don't believe anything in the first 3 radial bins
	#w[:3,:,:,:] = 0

	# go ahead and remove the table from memory
	del(table, norm)

	print('Number of knots used: ',[len(a) for a in knots])
	print("Beginning spline fit for timing table...")
	spline = glam.fit(z,w,centers,knots,order,smooth,penalties=penalties,monodim=3)

	print("Saving table to %s..." % prob_outputfile)
	spline.knots = [spline.knots[i] * axis_scale[i] for i
			    in range(0, len(spline.knots))]
	splinefitstable.write(spline, prob_outputfile)

	# clean up
	del(w,z,order,penalties,knots,spline)


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

