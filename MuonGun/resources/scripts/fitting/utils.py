
import tables, numpy, dashi

def histload(hdf, where):
	"""
	dashi.histload with optional normalization
	"""
	h = dashi.histload(hdf, where)
	node = hdf.getNode(where)
	if 'count' in node._v_attrs:
		h /= node._v_attrs['count']
		print('norm: %.1f' % node._v_attrs['count'])
	return h

# Monkey-patch some handy features into dashi
def points(self, differential=False):
	sp = dashi.scatterpoints.points2d()
	sp.x = self.bincenters
	sp.xerr = self.binwidths/2.
	sp.y = self.bincontent.copy()
	sp.yerr = self.binerror
	if differential:
		sp.y /= self.binwidths
		sp.yerr /= self.binwidths
	return sp
dashi.histogram.hist1d.points = points

from icecube.photospline import numpy_extensions # meshgrid_nd
def select_fit_bins(self, range, differential):
	if range is None:
		x = numpy_extensions.meshgrid_nd(*self._h_bincenters, lex_order=True)
		y = self.bincontent.copy() # copy to allow in-place differentiation later
		error = numpy.sqrt(self.squaredweights)
		slices = [slice(None)]*self.ndim
	else:
		if not len(range) == self.ndim:
			raise ValueError("range must be a sequence of length ndim")
		slices = []
		for r, centers in zip(range, self._h_bincenters):
			if isinstance(r, slice):
				slices.append(r)
			else:
				if not len(r) == 2:
					raise ValueError("range for each dimension must be a 2-tuple")
				idxs = numpy.searchsorted(centers, r)
				slices.append(slice(idxs[0], idxs[1]))
		centers = [c[sl] for c, sl in zip(self._h_bincenters, slices)]
		x = numpy_extensions.meshgrid_nd(*centers, lex_order=True)
		y = self.bincontent[slices].copy()
		error = numpy.sqrt(self.squaredweights[slices])
 
	for dim in differential:
		if dim < 0 or dim >= self.ndim:
			raise ValueError("This %d-d histogram doesn't have a dimension %d" % (self.ndim, dim))
		edges = self._h_binedges[dim][self._h_visiblerange[dim]][slices[dim]]
		shape = [1]*self.ndim
		shape[dim] = edges.size-1
		dx = numpy.diff(edges).reshape(shape)
		y /= dx
		error /= dx
	return x, y, error

dashi.histogram.histogram._select_fit_bins = select_fit_bins

def multi_leastsq(self, model, range=None, differential=[], **kwargs):
	"""
	The n-dimensional generalization of hist1d.leastsq

	:param range: A list of slices or (min, max) pairs giving the selection of
	              bins to include in the fit. There must be one entry for each
	              dimension in the histogram.
	:param differential: A list containing the indices of the dimensions that
	                     should be differentiated before fitting.
	"""
	x, y, error = self._select_fit_bins(range, differential)
	return dashi.fitting.leastsq(x, y, model, error, **kwargs)

dashi.histogram.histogram.leastsq = multi_leastsq

def multi_empty_like(self):
	return dashi.histogram.histogram(self.ndim, [e.copy() for e in self._h_binedges], self.labels, self.title)
dashi.histogram.histogram.empty_like = multi_empty_like

def load_group(fname, group='energy'):
	import os
	dirname, groupname = os.path.split(fname)
	if os.path.isfile(dirname):
		fname = dirname
		root = '/'+groupname
	else:
		root = ''
	with tables.openFile(fname) as hdf:
		import operator
		try:
			# For real fluxes, elements are stored individually 
			hdf.getNode('%s/PPlus/%s' % (root, group))
			elements = 'PPlus', 'He4Nucleus', 'Fe56Nucleus', 'N14Nucleus', 'Al27Nucleus'
			h = histload(hdf, '%s/%s/%s' % (root, elements[0], group))
			for e in elements[1:]:
			    h += histload(hdf, '%s/%s/%s' % (root, e, group))
		except tables.NoSuchNodeError:
			# For pseudo-fluxes, there is only one group
			h = histload(hdf, root+'/'+group)
	return h

def load_espec(fname, single=True, bias=50, transform=True):
	he = load_group(fname, 'energy')

	if single:
		# select single muons and sum over all radii
		he = he[:,:,1,:,:].project([0,1,3])
		eaxis = 2
	else:
		eaxis = 4
	
	# normalize w.r.t. E
	norm = he._h_bincontent.sum(axis=eaxis).reshape(he._h_bincontent.shape[:-1] + (1,))
	if not single:
		# normalize w.r.t. r
		norm = norm.sum(axis=eaxis-1).reshape(norm.shape[:-2] + (1,1))
	he._h_bincontent /= norm
	he._h_squaredweights /= norm*norm
	# convert to differential
	shape = (1,)*(eaxis) + (he.bincontent.shape[eaxis],)
	norm = numpy.diff(he._h_binedges[-1][1:-1]).reshape(shape)
	# also convert to differential in r
	if not single:
		shape = [1]*he.ndim
		shape[eaxis-1] = he.bincontent.shape[eaxis-1]
		shape = tuple(shape)
		norm = norm.repeat(shape[eaxis-1], axis=eaxis-1)
		if transform:
			norm *= numpy.diff(he._h_binedges[eaxis-1][1:-1]**2).reshape(shape)
		else:
			norm *= numpy.diff(he._h_binedges[eaxis-1][1:-1]).reshape(shape)
	vrange = he._h_visiblerange
	he._h_bincontent[vrange] /= norm
	he._h_squaredweights[vrange] /= norm*norm
	
	if transform:
		# convert energies to log-log
		he._h_bincontent[:] = numpy.log(he._h_bincontent) + bias
		he._h_binedges[-1][1:-1] = numpy.log(he._h_binedges[-1][1:-1])
	
		# convert zenith angles to cos(zenith), reversing the angular axis in the process
		he._h_binedges[0][1:-1] = numpy.cos(he._h_binedges[0][1:-1][::-1])
		# reverse through a temporary to avoid overwriting bits we want to read later
		rev = he._h_bincontent[::-1,:,:].copy()
		he._h_bincontent[:] = rev
		rev = he._h_squaredweights[::-1,:,:].copy()
		he._h_squaredweights[:] = rev
	
		# zero out non-finite weights
		mask = numpy.logical_not(numpy.isfinite(he._h_bincontent) & numpy.isfinite(he._h_squaredweights))
		he._h_bincontent[mask] = 0
		he._h_squaredweights[mask] = 0
	
	return he

def load_radial_distribution(fname, bias=50, transform=True):
	h = load_group(fname, 'radius')
	#normalize
	norm = h._h_bincontent.sum(axis=3).reshape(h._h_bincontent.shape[:-1] + (1,))
	h._h_bincontent /= norm
	h._h_squaredweights /= norm*norm

	# convert to differential
	shape = (1,)*(h.bincontent.ndim-1) + (h.bincontent.shape[-1],)
	if transform:
		norm = numpy.diff(h._h_binedges[-1][1:-1]**2).reshape(shape)
	else:
		norm = numpy.diff(h._h_binedges[-1][1:-1]).reshape(shape)
	vrange = h._h_visiblerange
	h._h_bincontent[vrange] /= norm
	h._h_squaredweights[vrange] /= norm*norm
	
	if transform:
	
		# parameterize log(dP/dR^2) as a function of R
		h._h_bincontent[:] = numpy.log(h._h_bincontent) + bias
	
		# convert zenith angles to cos(zenith), reversing the angular axis in the process
		h._h_binedges[0][1:-1] = numpy.cos(h._h_binedges[0][1:-1][::-1])
		# reverse through a temporary to avoid overwriting bits we want to read later
		rev = h._h_bincontent[::-1,:,:].copy()
		h._h_bincontent[:] = rev
		rev = h._h_squaredweights[::-1,:,:].copy()
		h._h_squaredweights[:] = rev
	
		# zero out non-finite weights
		mask = numpy.logical_not(numpy.isfinite(h._h_bincontent) & numpy.isfinite(h._h_squaredweights))
		h._h_bincontent[mask] = 0
		h._h_squaredweights[mask] = 0
	
	return h

def load_flux(fname, transform=False, bias=50):
	h = load_group(fname, 'multiplicity')
	norm = numpy.diff(numpy.cos(h._h_binedges[0][::-1])).reshape((h._h_bincontent.shape[0],) + (1,)*(h.ndim-1))	
	h._h_bincontent /= norm
	h._h_squaredweights /= norm*norm
	
	if h._h_binedges[2][1] > 0.5:
		# center bins on integers
		h._h_binedges[2] -= 0.5
	
		# iron out a stupid precision issue in the binning
		for sl in (slice(27,29), slice(54,56), slice(59,61)):
			h._h_bincontent[:,:,sl] = 0
			h._h_squaredweights[:,:,sl] = 0
	
	if transform:
		
		h._h_bincontent[h._h_bincontent < 1e-16] = 0
		
		# parameterize log(flux)
		h._h_bincontent[:] = numpy.log(h._h_bincontent) + bias
		
		# convert zenith angles to cos(zenith), reversing the angular axis in the process
		h._h_binedges[0][1:-1] = numpy.cos(h._h_binedges[0][1:-1][::-1])
		# reverse through a temporary to avoid overwriting bits we want to read later
		rev = h._h_bincontent[::-1,:,:].copy()
		h._h_bincontent[:] = rev
		rev = h._h_squaredweights[::-1,:,:].copy()
		h._h_squaredweights[:] = rev
	
		# zero out non-finite weights
		mask = numpy.logical_not(numpy.isfinite(h._h_bincontent) & numpy.isfinite(h._h_squaredweights))
		h._h_bincontent[mask] = 0
		h._h_squaredweights[mask] = 0
	
	return h

def pad_knots(knots, order=2):
	"""
	Pad knots out for full support at the boundaries
	"""
	pre = knots[0] - (knots[1]-knots[0])*numpy.arange(order, 0, -1)
	post = knots[-1] + (knots[-1]-knots[-2])*numpy.arange(1, order+1)
	return numpy.concatenate((pre, knots, post))

def colorize(sequence, cmap='jet'):
	from matplotlib.colors import Normalize
	from matplotlib import cm
	norm = Normalize(vmin=0, vmax=len(sequence)-1)
	cmap = getattr(cm, cmap)
	for i, s in enumerate(sequence):
		yield cmap(norm(i)), s

def plot_energy_slice(h, spline, slice_=(3,10,1,1), **kwargs):
	import pylab
	import dashi; dashi.visual()
	from icecube.photospline.glam.glam import grideval
	
	idx = slice_ + (slice(None),)
	sub = h[idx]
	idx = tuple([s-1 for s in slice_]) + (slice(None),)
	
	coords = []
	for i, s in enumerate(idx):
		axis = h._h_bincenters[i][s]
		try:
			len(axis)
			#axis = numpy.linspace(axis[0], axis[-1], 101)
			axis = numpy.linspace(0, 20, 101)
			# print axis
		except TypeError:
			axis = [axis]
		coords.append(axis)
		deg = (180*numpy.arccos(coords[0][0])/numpy.pi)
	sub.scatter(**kwargs)
	pylab.plot(coords[-1], grideval(spline, coords).flatten(), **kwargs)

def plot_radial_slice(h, spline, slice_=(3,10,1), **kwargs):
	import pylab
	import dashi; dashi.visual()
	from icecube.photospline.glam.glam import grideval
	
	idx = slice_ + (slice(None),)
	sub = h[idx]
	idx = tuple([s-1 for s in slice_]) + (slice(None),)
	
	coords = []
	for i, s in enumerate(idx):
		axis = h._h_bincenters[i][s]
		try:
			len(axis)
			#axis = numpy.linspace(axis[0], axis[-1], 101)
			axis = numpy.linspace(0, 250, 1001)
			# print axis
		except TypeError:
			axis = [axis]
		coords.append(axis)
	sub.scatter(**kwargs)
	print(grideval(spline, coords).flatten())
	pylab.plot(coords[-1], grideval(spline, coords).flatten(), **kwargs)
