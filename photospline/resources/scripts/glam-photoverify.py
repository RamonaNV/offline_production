from icecube.photospline.glam import glam, bspline
from icecube.photospline import splinetable, splinefitstable
from icecube.photospline.utils import TableSlice
from icecube.photospline.photonics import *

from optparse import OptionParser

from icecube.photospline.numpy_extensions import *

import sys, os
import numpy
import Gnuplot

usage = "usage: %prog [options] table.pt table.fits [table2.fits [ ... ] ]"
optparser = OptionParser(usage=usage)
optparser.add_option("-0", "--dim0", dest="dim0", type="string",
             help="How to plot histogram dimension 0 [x|y|i|<bin>]")
optparser.add_option("-1", "--dim1", dest="dim1", type="string",
             help="How to plot histogram dimension 1 [x|y|i|<bin>]")
optparser.add_option("-2", "--dim2", dest="dim2", type="string",
             help="How to plot histogram dimension 2 [x|y|i|<bin>]")
optparser.add_option("-3", "--dim3", dest="dim3", type="string",
             help='''How to plot histogram dimension 3 [x|y|i|<bin>]                               
                     x    : use this dimension as x variable in plot                               
                     y    : use this dimension as y variable in plot                               
                     i    : iterate over this dimension                                            
                     <bin>: slice at this bin number''')
optparser.add_option("--start",dest="start",type="int",default=0,
        help='Start iterating at this bin number')
optparser.add_option("--density",dest="density",type="int",default=1,
        help="Points to plot per table bin")
optparser.add_option("--dump",dest="dump",action="store_true",default=False,
        help="Iterate over dimension, dumping each plot to an EPS file.")
(opts, args) = optparser.parse_args()

if len(args) < 2:
    print(usage)
    sys.exit(1)

# Load original table. If it fails, try another format.
try:
	table = Table(args[0], mmap=False)
	geo = table.ph_header.geo
except:
	table = FITSTable.load(args[0])
	geo = table.header['geometry']
	
table.remove_nans_and_infinites()
table.normalize()

print("Loaded histogram with dimensions ",table.shape)

three_d = False
xdim = None
ydim = None
idim = None
axes = [opts.dim0, opts.dim1, opts.dim2, opts.dim3]
free_axes = list(range(table.ndim))

if 'y' in axes:
    ydim = axes.index('y') 
    if ydim >= table.ndim:
        print(ydim,"-> y: Table only has dimensions", list(range(table.ndim)))
        sys.exit()
    free_axes.remove(ydim)
    three_d = True

if 'x' in axes:
    xdim = axes.index('x')
    if xdim >= table.ndim:
        print(xdim,"-> x: Table only has dimensions", list(range(table.ndim)))
        sys.exit()
else:
    xdim = max(free_axes)
free_axes.remove(xdim)

if 'i' in axes:
    idim = axes.index('i')
    if idim >= table.ndim:
        print(idim,"-> i: Table only has dimensions", list(range(table.ndim)))
        sys.exit()
else:
    idim = max(free_axes)
free_axes.remove(idim)

class AxisDesc(object):
	def __init__(self, abbrev, desc, unit=None):
		self._abbrev = abbrev
		self._desc = desc
		self._unit = unit	
	@property
	def unit(self):
		if self._unit:
			return " [%s]" % self._unit
		else:
			return ""
	@property
	def description(self):
		return self._desc + self.unit
	def format(self, val):
		return "%s = %.2f%s" % (self._abbrev, val, self.unit)

axis_labels = [
	AxisDesc('rho', 'Perpendicular distance', 'm'),
	AxisDesc('theta', 'Observervation angle (wrt vertical)', 'degree'),
	AxisDesc('z', 'Parallel distance', 'm'),
	AxisDesc('t', 'Time delay', 'ns'),
]

# Label axes appropriately
if geo == Geometry.CYLINDRICAL:
	pass
elif geo == Geometry.SPHERICAL:
	axis_labels[0] = AxisDesc('r', 'Source-observer distance', 'm')
	axis_labels[2] = AxisDesc('cos(z)', 'Observervation angle (wrt source axis)', None)
else:
	print('Unknown table geometry %d!' % geo)

#print 'x:', xdim, '   y:', ydim, '   i:', idim, 'free:', free_axes

# Load the spline fit

splines = []
fit_labels = []
for file in args[1:]:
    spline = splinefitstable.read(file)
    splines.append(spline)
    fit_labels.append(os.path.splitext(file)[0])

# ---------------------------------------------------------------------------

def printdiffstats(a,b):
    print("Fit Statistics:")
    print("\tMaximum Deviation from Data:",numpy.max(numpy.abs(a - b)))
    print("\tRMS Deviation from Data:",numpy.sqrt(numpy.mean((a - b)**2)))
    print("\tMax Fractional Deviation from Data:",numpy.max(numpy.abs((a - b)/b)))
    print("\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs((a - b)/b)))


gp = Gnuplot.Gnuplot()

timing = len(spline.knots) > 3

for i,icenter in enumerate(table.bin_centers[idim]):
    if i < opts.start:
        continue
    title = "Slice at %s" % (axis_labels[idim].format(icenter))
    try:
        slices = [slice(None)]*4
        slices[idim] = i
        for d in free_axes:
            if axes[d] is not None:
                bin = int(axes[d])
            else:
                bin = int(len(table.bin_centers[d])/2.)
            slices[d] = bin
            title += ", %s" % (axis_labels[d].format(table.bin_centers[d][bin]))
        print("Building data set...")
        if len(splines) == 1:
            sample = TableSlice(table, spline, slices, opts.density).flatten()
        else:
            samples = [TableSlice(table, spline, slices, opts.density).flatten() for spline in splines]
            datacols = 2 if timing else 1
            stacks = tuple([samples[0]] + [sample[:,-datacols:] for sample in samples[1:]])
            sample = numpy.column_stack(stacks)

        gp.xlabel(axis_labels[xdim].description)
        gp.title(title)
        plots = []

        ndim = table.ndim
        xshift = numpy.zeros(len(sample))
        if timing:
            offset = 2
            # shift the raw CDF points to the right-hand bin edges
            xshift[::opts.density] = table.bin_widths[xdim]/2.0
        else:
            offset = 1
        if three_d:
            pass
            for plot in range(table.ndim+2,sample.shape[1],2):
                plots.append(Gnuplot.Data(sample[:,xdim],
                                          sample[:,ydim],
                                          sample[:,plot],
                                          title="%s CDF" % fit_labels[(plot - ndim - offset)/offset]))
            plots.append(Gnuplot.Data(sample[:,xdim], sample[:,ydim], sample[:,table.ndim],   title="Raw CDF"))
            #plots.append(Gnuplot.Data(sample[:,xdim], sample[:,ydim], sample[:,table.ndim+1],   title="Raw PDF"))
            gp.ylabel(axis_labels[ydim].description)
            gp.zlabel("CDF")
            gp.splot(*plots)
        else:
            showpdf = timing and xdim == 3
            for cidx, plot in enumerate(range(ndim+offset,sample.shape[1],offset)):
                if showpdf:
                    ccolor = cidx*2 + 1
                    pcolor = ccolor + 1
                else:
                    ccolor = cidx + 1
                plots.append(Gnuplot.Data(sample[:,xdim],
                                          sample[:,plot],
                                          title="%s CDF" % fit_labels[(plot - ndim - offset)/offset], axes='x1y1', with_='lines lc %d' % ccolor))
                if showpdf:
                    plots.append(Gnuplot.Data(sample[:,xdim],
                                          sample[:,plot+1],
                                          title="%s PDF" % fit_labels[(plot - ndim - offset)/offset], axes='x1y2', with_='lines lc %d' % pcolor))

            # mask out non-finite values to keep Gnuplot happy
            dmask = numpy.isfinite(sample[:,ndim])
            plots.append(Gnuplot.Data(sample[dmask,xdim]+xshift[dmask], sample[dmask,ndim],   title="Raw CDF", axes='x1y1', with_='points lc 1', every=opts.density))
            if showpdf:
                plots.append(Gnuplot.Data(sample[:,xdim], sample[:,ndim+1], title="Raw PDF", axes='x1y2', with_='points lc 3', every=opts.density))
            if timing:
                ylabel = "CDF"
                gp.ylabel(ylabel)
                ylabel = "PDF"
                gp("set y2label '%s'"%ylabel)
                #gp("set y2tics 0, 0.02")
                #gp("set y2range [*:*]")
                #gp("set ytics nomirror")
            else:
                gp.ylabel('log Amplitude [PE/m track]')
            
            yvals = sample[:,ndim]
            ymin = yvals[numpy.isfinite(yvals)].min()
            if ymin == 0:
                ymin += 1e-20
            ymax = yvals.max()
            gp('set style line 9 lt rgb "#D0D0D0" lw 0.1')
            gp('set style arrow 2 nohead ls 9')
            for k in spline.knots[xdim]:
                gp('set arrow from %e,%e to %e,%e as 2'%(k,ymin,k,ymax))
            gp.plot(*plots)

        printdiffstats(sample[::opts.density,ndim], sample[::opts.density,ndim+offset])
        if opts.dump:
            if not three_d:
                gp.set(yrange=(0,1.2))
                gp.set_boolean('grid',True)
            gp.hardcopy('cdf_slice_%s_%.3d_%s_%f.eps'%(axis_vars[idim], i, axis_vars[d], table.bin_centers[d][bin]),eps=True,color=True)
        else:
            gp.interact()
            #raw_input('Press Enter to continue')
    except AssertionError:
        print('Skipping   slice at  (no raw data)' % i)
