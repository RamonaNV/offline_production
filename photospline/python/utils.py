import copy
import numpy
from . import splinefitstable
from .glam import glam, bspline

def pad_knots(knots, order=2):
	"""
	Pad knots out for full support at the boundaries
	"""
	pre = knots[0] - (knots[1]-knots[0])*numpy.arange(order, 0, -1)
	post = knots[-1] + (knots[-1]-knots[-2])*numpy.arange(1, order+1)
	return numpy.concatenate((pre, knots, post))

class TableSlice(object):
	"""A slice of a photonics table, with spline CDF and PDF evaluates."""
	table_pdf = None
	table_cdf = None
	spline_pdf = None
	spline_cdf = None
	edges = None
	centers = None
	def __init__(self, table, spline, slices, density = 1):
		self.table = table
		self.spline = spline
		self.slices = slices[:table.values.ndim]
		self.density = density
		self.make_grid()
		if spline.level == 2:
			if len(spline.knots) == 3:
				norm = True
				is_log = False
			else:
				norm = False
				is_log = True
		else:
			if len(spline.knots) == 4:
				norm = True
				is_log = False
			else:
				norm = False
				is_log = True
		self.slice(is_log, norm)
		self.eval(is_log)

	def collapse(self):
		"""Collapse non-extended dimensions, returning 1 or 2-d arrays
		   for use with matplotlib.
		"""
		new = copy.copy(self)
		ext_dims = [i for i in range(len(self.centers)) if len(self.centers[i]) > 1]
		new.edges = [self.edges[i] for i in ext_dims]
		new.centers = [self.centers[i] for i in ext_dims]
		return new
		
	def flatten(self):
		"""Flatten grid to columns for use with Gnuplot"""
		coords = numpy.meshgrid_nd(*self.centers,**{'lex_order':True})
		coords = numpy.column_stack([arr.reshape(arr.size) for arr in coords])
		bigdat = numpy.column_stack((coords,self.table_cdf.flatten()))
		if self.table_pdf is not None:
			bigdat = numpy.column_stack((bigdat,self.table_pdf.flatten()))
		bigdat = numpy.column_stack((bigdat,self.spline_cdf.flatten()))
		if self.spline_pdf is not None:
			bigdat = numpy.column_stack((bigdat,self.spline_pdf.flatten()))
		return bigdat

	def make_grid(self, density = 1):
		slices = self.slices
		centers = [bins[_slice] for bins,_slice in zip(self.table.bin_centers,slices)]
		for i,sl in enumerate(slices):
			if isinstance(sl,int):
				centers[i] = numpy.array([centers[i]])
			elif self.density > 1: # add extra grid points in between the bin centers
				gaps = numpy.diff(centers[i])
				extras = []
				for j in range(1,self.density):
					scale = j/float(self.density)
					extras.append(centers[i][:-1]+scale*gaps)
				centers[i] = numpy.concatenate(tuple([centers[i]] + extras))
				centers[i].sort()
		self.centers = centers

		widths = [bins[_slice] for bins,_slice in zip(self.table.bin_widths,slices)]
		for i,sl in enumerate(slices):
			if isinstance(sl,int):
				widths[i] = numpy.array([widths[i]])
			elif self.density > 1:
				# subdividing the widths is much easier!
				rep = self.density*numpy.ones(widths[i].size, dtype=int)
				rep[-1] = 1
				widths[i] = (widths[i]/self.density).repeat(rep)
		edges = [c - w/2.0 for c,w in zip(centers,widths)]
		edges = [numpy.append(e, c[-1]+w[-1]/2.0) for e,c,w in zip(edges,centers,widths)]
		self.edges = edges

		if len(self.spline.knots) == 4:
			# XXX: evaluate CDF at right edge of the time bin
			self.centers[3] = self.edges[3][1:]

	def slice(self, is_log = False, norm = True):
		timing = len(self.spline.knots) == 4
		table_pdf = None
		table_cdf = None
		slices = self.slices
		table = self.table
		if self.table.normed and self.slices[-1] == slice(None): # already normed
			table_cdf = self.table.values[slices].copy()
			if timing:
				print(table_cdf.shape, table.bin_widths[-1].size)
				table_pdf = numpy.append([0],
				    numpy.diff(table_cdf, axis=-1))/table.bin_widths[-1]
		elif self.table.normed:
			# already normed, but we're slicing perpendicular to time.
			# skip the PDF.
			table_cdf = self.table.values[slices].copy()
		elif len(self.spline.knots) == 3:
			# abs spline, just sum amplitudes
			tslices = list(self.slices)[:3]
			tslices += [slice(None)]
			bigslice = self.table.values[tslices]
			table_cdf = numpy.sum(bigslice, axis=-1)
		elif timing and self.slices[-1] == slice(None): # we're slicing in time, compute CDF for slice
			table_slice = self.table.values[slices]
			#table_cdf = numpy.cumsum(table_slice*table.bin_widths[3], axis=-1)
			table_cdf = numpy.cumsum(table_slice, axis=-1)
			#table_pdf = table_slice
			table_pdf = table_slice/table.bin_widths[3]
			if norm:
				normslices = [slice(None)]*len(table_cdf.shape)
				normslices[-1] = -1
				newshape = list(table_cdf.shape)
				newshape[-1] = 1
				normval = table_cdf[tuple(normslices)].reshape(tuple(newshape))
				table_cdf /= normval
				table_pdf /= normval
		else:
			# we're slicing perpendicular to time, so slice in the dim-t plane
			# , compute the CDF, then slice in t
			tslices = list(self.slices)
			if timing:
				tslices[-1] = slice(None)
			else:
				tslices += [slice(None)]
			print(tslices)
			bigslice = self.table.values[tslices]
			table_cdf = numpy.cumsum(bigslice, axis=-1)
			# remove integer indexes, since those dimensions are already gone
			tslices = [t for t in tslices if not isinstance(t,int)]
			# now slice at the time we want
			tslices[-1] = slices[-1]
			nslice = [slice(None)]*table_cdf.ndim
			nslice[-1] = -1
			normval = table_cdf[nslice]
			table_cdf = table_cdf[tslices]
			if timing:
				# careful! table_pdf is a view into table.values
				table_pdf = table.values[slices].copy()
			else:
				table_cdf = normval
			if norm:
				table_cdf /= normval
				table_pdf /= normval

		if self.density > 1: # insert zeros into table values
			expanded_shape = tuple([s + (self.density-1)*(s-1) for s in table_cdf.shape])
			insert_slices = [slice(None,None,self.density)]*len(table_cdf.shape)
			t_cdf = numpy.zeros(expanded_shape)
			t_cdf[insert_slices] = table_cdf
			table_cdf = t_cdf
			if timing:
				t_pdf = numpy.zeros(expanded_shape)
				t_pdf[insert_slices] = table_pdf
				table_pdf = t_pdf
		self.table_cdf = table_cdf
		self.table_pdf = table_pdf
		
	def eval_cdf(self, is_log):
		vals = glam.grideval(self.spline, self.centers)
		shape = vals.shape
		vals = vals[numpy.isfinite(vals)] + self.spline.bias
		shape = tuple([shape[i] for i in range(len(shape)) if shape[i] > 1])
		vals = vals.reshape(shape)
		if is_log:
			vals = numpy.exp(vals)
		self.spline_cdf = vals
	def eval_pdf(self, is_log):
		
		splfuncs = [bspline.bspline]*4
		splfuncs[-1] = bspline.bspline_deriv
		deriv_basis = [bspline.splinebasis(self.spline.knots[i], self.spline.order[i],self.centers[i],
                    self.spline.periods[i],spline=splfuncs[i]) for i in range(0,len(self.spline.knots))]
		pdf_vals = glam.grideval(self.spline, self.centers, bases = deriv_basis)
		shape = pdf_vals.shape
		pdf_vals = pdf_vals[numpy.isfinite(pdf_vals)] #+ spline.bias
		if is_log:
			# chain rule!
			pdf_vals *= self.spline_cdf
		shape = tuple([shape[i] for i in range(len(shape)) if shape[i] > 1])
		pdf_vals = pdf_vals.reshape(shape)
		self.spline_pdf = pdf_vals
		
		
	def eval(self, is_log = False):
		"""is_log => fit is in log-space"""
		self.eval_cdf(is_log)
		# only calculate PDFs for timing fits
		if (len(self.spline.knots)) == 4:
			self.eval_pdf(is_log)

