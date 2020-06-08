
import numpy, warnings, struct

class MetaHead(object):
	"""
	Parse a packed representation of MetaHead_type from photonics.h.
	"""
	_struct = struct.Struct('<5s8B32s40s')
	size = _struct.size
	def __init__(self, buf):
		v = self._struct.unpack(buf)
		stringify = lambda s: s[:s.index('\0')]
		self.File_Format_Label = stringify(v[0])
		self.File_Format_Revision = v[1]
		self.bitsys = v[3]
		self.endian = v[4]
		self.level  = v[5]
		self.Photonics_Version = stringify(v[9])
		
	def pack(self):
		v = (self.File_Format_Label, self.File_Format_Revision, 0, self.bitsys,
		    self.endian, self.level, 0, 0, 0, self.Photonics_Version, "")
		return self._struct.pack(*v)
		
class Header(object):
	"""
	Parse a packed representation of Header_type from photonics.h.
	"""
	_struct = struct.Struct('<100s6f3i6f7i25f2i1f1i2q2i')
	size = _struct.size
	def __init__(self, fh):
		v = self._struct.unpack(fh.read(self.size))
		self.MetaHead = MetaHead(v[0][:-15])
		self.sphere_frac = v[1]
		self.impsampl_Lf = v[2]
		self.impsampl_Tf = v[3]
		self.ref_np      = v[5]
		self.ref_ng      = v[6]
		self.record_errors = bool(v[7])
		self.source_type = v[8]
		self.extended_source = bool(v[9])
		self.step        = v[10]
		self.e           = v[11]
		self.volume      = v[12]
		self.angle       = v[13]
		self.source_rz   = (v[14], v[15])
		self.geo         = v[16]
		self.n           = v[17:23]
		self.limits = []
		self.maxes = []
		pos = 23
		for i in range(6):
			self.limits.append(v[pos:pos+2])
			pos += 2
		for i in range(6):
			self.maxes.append(v[pos:pos+2])
			pos += 2
		self.depth       = v[47]
		self.d_scale     = v[48]
		self.t_scale     = v[49]
		self.lambda_     = v[50]
		self.efficiency  = v[51]
		self.n_photon    = v[52]
		self.n_entries   = v[53]
		self.refraction_mode = v[54]
		
	def write(self, fh):
		v = [self.MetaHead.pack(), self.sphere_frac, self.impsampl_Lf, self.impsampl_Tf, 0,
		    self.ref_np, self.ref_ng, self.record_errors, self.source_type, self.extended_source,
		    self.step, self.e, self.volume, self.angle]
		v += list(self.source_rz)
		v += [self.geo]
		v += list(self.n)
		for l in self.limits:
			v += list(l)
		for l in self.maxes:
			v += list(l)
		v += [self.depth, self.d_scale, self.t_scale, self.lambda_, self.efficiency,
		    self.n_photon, self.n_entries, self.refraction_mode, 0]
		fh.write(self._struct.pack(*v))

assert(Header.size == 328)
		
class Efficiency(object):
	"""Normalization types from photonics.h"""
	NONE         = 0x00
	RECEIVER     = 0x01 
	SOURCE       = 0x02
	WAVELENGTH   = 0x04
	AREA         = 0x08
	VOLUME       = 0x10
	N_PHOTON     = 0x20
	EMISSION     = 0x40
	USER_DEFINED = 0x80
	DIFFERENTIAL = 0x100
	N_EVENT      = 0x200

class Geometry(object):
	SPHERICAL   = 1
	CYLINDRICAL = 2
	CUBIC       = 3

class GeometryType(object):
	POINTSOURCE = 0
	INFINITEMUON = 1

class Parity(object):
	EVEN	    = 0
	ODD	    = 1

# Class for reading IceCube photonics tables
class Table(object):
	level = -1

	table       = None
	weights     = None
	bin_centers = None
	bin_widths  = None

	is_integral = False

	filename    = None

	# Constructor. Creates instance and optionally opens pt file.
	def __init__(self, filename=None, normalize=True, mmap=True):
		if filename is not None:
			self.open_file(filename, mmap=mmap)
		if normalize:
			try:
				self.normalize()
			except:
				pass

	@classmethod
	def stack(cls, outfile, *fnames):
		import shutil
		shutil.copy(fnames[0], outfile)
		target = cls()
		target._read_standalone(outfile, mmap=True, mode='r+')
		for fn in fnames[1:]:
			piece = cls()
			piece._read_standalone(fn, mmap=True, mode='r')
			target.values += piece.values
			if piece.ph_header.record_errors and target.ph_header.record_errors:
				target.weights += piece.weights
			target.ph_header.n_photon += piece.ph_header.n_photon
		with open(outfile, 'r+') as fh:
			fh.seek(0)
			target.ph_header.write(fh)

	# Checks consistency of loaded tables.
	# For now, only checks shapes of various arrays.
	def table_shape_consistent(self):
		shapes = set()

		shapes.add(self.values.shape)
		if self.weights is not None:
			shapes.add(self.weights.shape)
		shapes.add(tuple([len(i) for i in self.bin_centers]))
		shapes.add(tuple([len(i) for i in self.bin_widths]))
		if len(shapes) > 1:
			return 0

		return 1

	# Normalize to absolute bin amplitudes
	def normalize(self):
		eff = self.header['efficiency']
		if not (eff & Efficiency.N_PHOTON):
			# This table still contains raw weights from photomc
			self.values /= self.header['n_photon']
			self.weights /= self.header['n_photon']
			eff = eff | Efficiency.N_PHOTON
		if (eff & Efficiency.DIFFERENTIAL):
			# Someone has made this a dP/dt table. Undo their work.
			if self.values.ndim != 4:
				raise ValueError("This table is weird, man.")
			shape = [1]*len(self.values.shape)
			shape[-1] = self.values.shape[-1]
			dt = self.bin_widths[-1].reshape(shape)
			self.values *= dt
			eff = eff & ~Efficiency.DIFFERENTIAL

		self.header['efficiency'] = eff

	# Returns number of dimensions in table.
	@property
	def ndim(self):
		return len(self.shape)

	# Returns shape of table.
	@property
	def shape(self):
		if not self.table_shape_consistent():
			raise Exception('Shape consistency check failed')

		return self.values.shape

	def remove_nans_and_infinites(self, dovalues=True, doweights=True):
		if self.weights is not None and doweights:
			self.weights[numpy.logical_not( \
			    numpy.isfinite(self.values))] = 0
		if dovalues:
			self.values [numpy.logical_not( \
			    numpy.isfinite(self.values))] = 0

	@property
	def normed(self):
		"""Has this table been normalized?"""
		if self.values.ndim == 4:
			normval = self.values[:,:,:,-1]
			if (normval[(normval > 0) & \
			    numpy.isfinite(normval)] == 1).all():
				return True
			else:
				return False
		else:
			return True
			
	def _read_standalone(self, filename, mmap, mode='r'):
		"""
		Read using standalone implementation.
		"""
		
		header = Header(open(filename))
		
		if header.MetaHead.level == 2:
			raise ValueError("I don't know how to read level-2 tables!")
		self.values = numpy.squeeze(numpy.memmap(filename, shape=header.n, dtype=numpy.float32, offset=Header.size, mode=mode))
		if header.record_errors:
			offset = Header.size + self.values.itemsize*self.values.size
			self.weights = numpy.squeeze(numpy.memmap(filename, shape=header.n, dtype=numpy.float32, offset=offset, mode=mode))
		else:
			self.weights = numpy.zeros(self.values.shape, dtype=numpy.float32)

		# In keeping with a convention established by photo2numpy,
		# tables are either mmap'd in single precision or read in
		# to memory completely in double precision
		if not mmap:
			self.values = self.values.astype(numpy.float64)
			self.weights = self.weights.astype(numpy.float64)
		else:
			warnings.warn("Memory-mapped tables are single-precision. You have been warned.");

		self.bin_centers = []
		self.bin_widths = []
		
		trafos = [lambda a: a]*len(header.limits)
		itrafos = [lambda a: a]*len(header.limits)
		if header.geo == Geometry.SPHERICAL:
			trafos[2] = lambda a: -numpy.cos(numpy.pi*a/180.)
		if header.d_scale == 2:
			trafos[0] = numpy.sqrt
			itrafos[0] = lambda a: a**2
		if header.t_scale == 2:
			trafos[-1] = numpy.sqrt
			itrafos[-1] = lambda a: a**2
		
		for i in range(len(header.limits)):
			steps = header.n[i]+1
			if steps == 2:
				continue
			lo, hi = list(map(trafos[i], header.limits[i]))
			edges = itrafos[i](numpy.linspace(lo, hi, steps))
			self.bin_centers.append(0.5*(edges[1:]+edges[:-1]))
			self.bin_widths.append(numpy.diff(edges))
		
		self.ph_header = header
		
		# Add compatibility 
		self.level = header.MetaHead.level	
		self.header = {
			'n_photon'  : header.n_photon,
			'efficiency': header.efficiency,
			'geometry'  : header.geo,
			'zenith'    : header.angle,
			'z'         : header.depth,
			'n_group'   : header.ref_ng,
			'n_phase'   : header.ref_np,
		}

	def open_file(self, filename, convert=False, mmap=False):
		
		self._read_standalone(filename, mmap)

		# Level 2 tables get an extra, random element here for
		# some reason
		if len(self.bin_centers) > self.values.ndim:
			self.bin_centers = self.bin_centers[0:self.values.ndim]
		if len(self.bin_widths ) > self.values.ndim:
			self.bin_widths  = self.bin_widths [0:self.values.ndim]
            
		# Check consistency of table shapes and derive type
		ndim = self.ndim
		if ndim == 3:
			self.is_integral = True;

		# Convert to standard format unless user doesn't want this
		if convert:
			self.convert_to_level1()

		return 1

	def convert_to_level1(self):
		if self.level == 0 or self.level == 1:
			return 1

		# For level 2, some axes are reversed.
		if self.level == 2:
			self.values = numpy.rollaxis(self.values, 0, 3)
			if self.weights is not None:
				self.weights = \
				    numpy.rollaxis(self.weights, 0, 3)
			self.bin_centers[2], self.bin_centers[0], \
			    self.bin_centers[1] = self.bin_centers[0], \
			    self.bin_centers[1], self.bin_centers[2]
			self.bin_widths[2],  self.bin_widths[0], \
			    self.bin_widths[1] = self.bin_widths[0], \
			    self.bin_widths[1],  self.bin_widths[2]
            
			from math import pi
			self.bin_centers[1][:] *= 180./pi
			self.bin_widths[1][:]  *= 180./pi
            
			self.level = 1

			return 1
            
		print("Don't know how to convert table with level", self.level)
		return 0
            
	def convert_to_level2(self):
		if self.level == 2:
			return 1

		# For level 0/1, some axes are reversed.
		if self.level == 0 or self.level == 1:
			self.values =       numpy.rollaxis(self.values, 2, 0)
			if self.weights is not None:
				self.weights = numpy.rollaxis(self.weights, \
				    2, 0)
			self.bin_centers[0], self.bin_centers[1], \
			    self.bin_centers[2] = self.bin_centers[2], \
			    self.bin_centers[0], self.bin_centers[1]
			self.bin_widths[0], self.bin_widths[1], \
			    self.bin_widths[2] = self.bin_widths[2], \
			    self.bin_widths[0], self.bin_widths[1]

			from math import pi
			self.bin_centers[1][:] *= pi/180.
			self.bin_widths[1][:]  *= pi/180.

			self.level = 2
            
			return 1
            
		print("Don't know how to convert table with level", self.level)
		return 0

	def mirror(self,n_rho=0,n_phi=0):
		"""Extend table to rho < 0 and 180 < phi < 360. This may be useful for surpressing edge effects while fitting."""

		if n_rho == 0 and n_phi == 0:
			return None

		if abs(self.bin_widths[1].sum() - 180) > 1e-12:
			raise ValueError("Only half-cylindrical tables can \
			    be mirrored. Perhaps mirror() has already been \
			    called?")

		## XXX only phi mirroring for now
		new_shape = list(self.values.shape)
		new_shape[0] += n_rho
		new_shape[1] += 2*n_phi

		target_slice = [slice(None)]*self.values.ndim
		source_slice = [slice(None)]*self.values.ndim
		target_slice[0] = slice(n_rho, None)
		target_slice[1] = slice(n_phi, -n_phi)

		# copy values into expanded array
		new_values = numpy.empty(new_shape)
		new_values[target_slice] = self.values

		# replace old values with expanded version
		del self.values
		self.values = new_values

		# copy weights into expanded array
		new_weights = numpy.empty(new_shape)
		new_weights[target_slice] = self.weights

		# replace weights
		del self.weights
		self.weights = new_weights

		# replace bin centers and widths
		for lst in (self.bin_centers, self.bin_widths):
			for i in (0,1):
				new = numpy.empty(new_shape[i])
				new[target_slice[i]] = lst[i]
				lst[i] = new

		# mirror left edge
		source_slice[1] = [2*n_phi - 1 - i for i in range(n_phi)]
		target_slice[0] = slice(None)
		target_slice[1] = list(range(n_phi))
		for array in (self.values, self.weights):
			array[target_slice] = array[source_slice]
		for lst in (self.bin_centers, self.bin_widths):
			lst[1][target_slice[1]] = -(lst[1][source_slice[1]])

		# mirror right edge
		source_slice[1] = [-(2*n_phi - i) for i in range(n_phi)]
		target_slice[1] = [-(i+1) for i in range(n_phi)]
		for array in (self.values, self.weights):
			array[target_slice] = array[source_slice]
		for lst in (self.bin_centers, self.bin_widths):
			lst[1][target_slice[1]] = 360 - lst[1][source_slice[1]]

		# mirror radial slices
		# negative radii are mirrored, so in reverse order
		source_slice[0] = list(range(2*n_rho - 1, n_rho - 1, -1))
		target_slice[0] = list(range(n_rho))
		for lst in (self.bin_centers, self.bin_widths):
			lst[0][target_slice[0]] = -(lst[0][source_slice[0]])

		# mirror the radial slice at each azimuth to negative radii
		for i in range(self.bin_centers[1].size):
			# find the opposite slice
			opposite = 180 + self.bin_centers[1][i]
			if opposite > 180: opposite -= 2*(180 - opposite)
			elif opposite < 0: opposite *= -1

			mcenter = abs(opposite - self.bin_centers[1]).argmin()
			source_slice[1] = mcenter
			target_slice[1] = i
			for array in (self.values, self.weights):
				array[target_slice] = array[source_slice]
	
		return None

class FITSTable(Table):
	"""
	Same content as a photonics table, but using the FITS-based file format
	produced by the clsim tabulator.
	"""
	
	# A default header. This contains the same keys as the one created in photo2numpy from photospline.
	empty_header = {
		'n_photons':         0,
		'efficiency':        Efficiency.NONE,
		'geometry':          Geometry.SPHERICAL,
		'parity':            Parity.EVEN,
		'zenith':            0.,
		'azimuth':           0.,
		'z':                 0.,
		'energy':            0.,
		'type':              0,
		'level':             1,
		'n_group':           numpy.nan,
		'n_phase':           numpy.nan,
	}
	
	def __init__(self, binedges, values, weights, header=empty_header):
		self.bin_edges = binedges
		self._visible_range = [slice(1,-1)]*len(binedges)
		shape = tuple((len(edges)+1 for edges in binedges))
		# add under and overflow bins if not present
		if values.shape == tuple((len(edges)-1 for edges in binedges)):
			full_values = numpy.zeros(shape)
			full_values[self._visible_range] = values
			values = full_values
			if weights is not None:
				full_weights = numpy.zeros(shape)
				full_weights[self._visible_range] = weights
				weights = full_weights
		assert values.shape == shape, "Data array has the correct shape"
		self._values = values
		self._weights = weights
		self.bin_centers = [(edges[1:]+edges[:-1])/2. for edges in self.bin_edges]
		self.bin_widths = [numpy.diff(edges) for edges in self.bin_edges]
		self.header = header
	
	@property
	def values(self):
		return self._values[self._visible_range]
	
	@property
	def weights(self):
		if self._weights is None:
			return None
		else:
			return self._weights[self._visible_range]
	
	def __getitem__(self, slice_):
		for i in slice_:
			if not (isinstance(i, int) or i == slice(None)):
				# print(slice_)
				print(i, isinstance(i, int), i == slice(None))
				raise ValueError("Only dimension-reducing slices are implemented")
		edges = [e for i, e in enumerate(self.bin_edges) if slice_[i] == slice(None)]
		if self.weights is None:
			w = None
		else:
			w = self.weights[slice_]
		return FITSTable(edges, self.values[slice_], w, self.header)
	
	def __iadd__(self, other):
		self.raise_if_incompatible(other)
		self._values += other._values
		if self._weights is not None:
			self._weights += other._weights
		self.header['n_photons'] += other.header['n_photons']
		if 'n_events' in self.header:
			self.header['n_events'] += other.header['n_events']
		
		return self

	def __idiv__(self, num):
		return self.__imul__(1./num)
		
	def __imul__(self, num):
		self._values *= num
		if self.weights is not None:
			self._weights *= num*num
		
		return self
		
	def raise_if_incompatible(self, other):
		"""
		Check for generalized brokenness.
		"""
		if not isinstance(other, self.__class__):
			raise TypeError("Can't combine a %s with this %s" % (other.__class__.__name__, self.__class__.__name__))
		if self._values.shape != other._values.shape:
			raise ValueError("Shape mismatch in data arrays!")
		nans = self._values.size - numpy.isfinite(self._values).sum()
		if nans != 0:
			raise ValueError("This table has %d NaN values. You might want to see to that.")
		nans = other._values.size - numpy.isfinite(other._values).sum()
		if nans != 0:
			raise ValueError("Other table has %d NaN values. You might want to see to that.")
		for k, v in self.header.items():
			if k in ('n_photons', 'n_events'):
				continue
			if other.header[k] != v:
				raise ValueError("Can't combine tables with %s=%s and %s" % (k, v, other.header[k]))
		
	def normalize(self, kind='photon'):
		"""
		Normalize the table. If *kind* is 'photon', normalize such that the
		entries in the table are number of PE detected per Cherenkov photon
		emitted between 300 and 600 nm, as in Photonics. If *kind* is 'event',
		normalize such that the entries in the table are the number of PE
		detected per event (e.g. per minimum-ionizing muon).
		"""
		if kind == 'photon':
			assert not self.header['efficiency'] & Efficiency.N_EVENT
			if not self.header['efficiency'] & Efficiency.N_PHOTON:
				self /= self.header['n_photons']
				self.header['efficiency'] |= Efficiency.N_PHOTON
		elif kind == 'event':
			assert not self.header['efficiency'] & Efficiency.N_PHOTON
			if not self.header['efficiency'] & Efficiency.N_EVENT:
				self /= self.header['n_events']
				self.header['efficiency'] |= Efficiency.N_EVENT
		else:
			raise ValueError("Unknown normalization type '%s'" % kind)
		
	def save(self, fname, overwrite=False):
		try:
			import pyfits
		except ImportError:
			import astropy.io.fits as pyfits
		import os
		
		if os.path.exists(fname):
			if overwrite:
				os.unlink(fname)
			else:
				raise IOError("File '%s' exists!" % fname)
		
		data = pyfits.PrimaryHDU(self._values)
		data.header.set('TYPE', 'Photon detection probability table')
		
		for k, v in self.header.items():
			# work around 8-char limit in FITS keywords
			tag = 'hierarch _i3_' + k
			data.header.set(tag, v)
		
		hdulist = pyfits.HDUList([data])
		
		if self._weights is not None:
			errors = pyfits.ImageHDU(self._weights, name='ERRORS')
			hdulist.append(errors)
		
		for i in range(self._values.ndim):
			edgehdu = pyfits.ImageHDU(self.bin_edges[i],name='EDGES%d' % i)
			hdulist.append(edgehdu)
			
		hdulist.writeto(fname)
		
	@classmethod
	def load(cls, fname):
		try:
			import pyfits
		except ImportError:
			import astropy.io.fits as pyfits
		
		hdulist = pyfits.open(fname)
		data = hdulist[0]
		values = data.data
		binedges = []
		for i in range(values.ndim):
			binedges.append(hdulist['EDGES%d' % i].data)
		
		try:
			weights = hdulist['ERRORS'].data
		except KeyError:
			weights = None
		
		header = dict()
		for k in map(str.lower, data.header.keys()):
			if k.startswith('_i3_'):
				header[k[4:]] = data.header[k]
			
		return cls(binedges, values, weights, header)

	@classmethod
	def stack(cls, outfile, *fnames):
		import os
		assert(not os.path.exists(outfile))
		target = cls.load(fnames[0])
		for fn in fnames[1:]:
			target += cls.load(fn)
		target.save(outfile)
		

def melonball(table, weights = None, radius = 1):
	"""Set weights inside a given radius to zero."""
	if table.header['geometry'] != Geometry.CYLINDRICAL:
		raise ValueError("Can't handle non-cylindrical tables")
	Rho, Z = numpy.meshgrid_nd(table.bin_centers[0], table.bin_centers[2], lex_order=True)
	mask = Rho**2 + Z**2 < radius**2
	if weights is None:
		weights = table.weights
	shape = weights.shape
	for i in range(shape[1]):
		if weights.ndim == 3:
			weights[:,i,:][mask] = 0
		else:
			for j in range(shape[3]):
				weights[:,i,:,j][mask] = 0
				
def scalp(table, weights = None, low = -820, high = 820):
	"""Set weights outside of the given depth range to zero."""
	
	geo = table.header['geometry']
	depth = table.header['z']
	zenith = table.header['zenith']*numpy.pi/180.0
	if geo == Geometry.CYLINDRICAL:
		Rho, Phi, L = numpy.meshgrid_nd(*table.bin_centers[:3], lex_order=True)
	elif geo == Geometry.SPHERICAL:
		R, Phi, CosPolar = numpy.meshgrid_nd(*table.bin_centers[:3], lex_order=True)
		L = R*CosPolar
		Rho = numpy.sqrt(R**2 - L**2)
		del R, CosPolar
	else:
		raise ValueError("Unknown geometry type %d" % geo)
	Phi *= (numpy.pi/180.0)
		
	z = L*numpy.cos(zenith) + Rho*numpy.sin(zenith)*(numpy.cos(Phi) + numpy.sin(Phi))
	mask = (z > low)&(z < high)
	
	del Rho, Phi, L, z
	
	if weights is None:
		weights = table.weights
	shape = weights.shape
	ddim = len(shape) - len(mask.shape)
	if ddim != 0:
		mask = mask.reshape(mask.shape + (1,)*ddim)
	
	weights[mask] = 0
	
