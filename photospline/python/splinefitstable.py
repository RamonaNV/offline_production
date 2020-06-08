from . import splinetable
try:
	import pyfits
except ImportError:
	import astropy.io.fits as pyfits
import numpy

def write(table,path):
	"""
	Write a SplineTable to disk as a FITS file.
	
	:param table: the SplineTable to be written
	:param path: where to save the FITS file.
	
	.. warning:: pyfits will fail to write the file if it already exists.
	"""
	data = pyfits.PrimaryHDU(table.coefficients)
	data.header.set('TYPE','Spline Coefficient Table')

	if getattr(table.order,'__iter__',False):
		for i in range(0,len(table.order)):
			data.header.set('ORDER%d' % i,table.order[i],
			    'B-Spline Order')
	else:
		data.header.set('ORDER',table.order,'B-Spline Order')

	for i in range(0,len(table.periods)):
		data.header.set('PERIOD%d' % i,table.periods[i])

	base_keys = set(['coefficients', 'order', 'knots', 'extents', 'periods'])
	for k in dir(table):
		if k.startswith('_') or k in base_keys:
			continue
		data.header.set(k.upper(), getattr(table, k))

	hdulist = pyfits.HDUList([data])

	for i in range(0,len(table.knots)):
		knothdu = pyfits.ImageHDU(table.knots[i],name='KNOTS%d' % i)
		hdulist.append(knothdu)
	
	extents = []
	for i in range(0,len(table.knots)):
		if getattr(table.order,'__iter__',False):
			order = table.order[i]
		else:
			order = table.order

		if i < len(table.extents):
			ext = list(table.extents[i])
		else:
			ext = [table.knots[i][order], table.knots[i][-order-1]]
		extents += ext
	extenthdu = pyfits.ImageHDU(numpy.array(extents, dtype=float), name='EXTENTS')
	hdulist.append(extenthdu)

	hdulist.writeto(path)

def read(path, memmap=False):
	"""
	Read a SplineTable from a FITS file on disk
	
	:param path: the filesystem path to read from
	:param memmap: memmap the underlying fits file this causes the underlying file handle to remain open indefinitely
	:returns: SplineTable - the spline surface stored in the given file
	"""
	
	file = pyfits.open(path, memmap=memmap)
	table = splinetable.SplineTable()

	data = file[0]
	table.coefficients = data.data
	table.periods = []
	table.knots = []
	for i in range(0,table.coefficients.ndim):
		try:
			table.periods.append(data.header['PERIOD%d' % i])
		except KeyError:
			table.periods.append(0.)
		table.knots.append(file['KNOTS%d' % i].data)


	try:
		table.order = data.header['ORDER']
		order = [table.order]*table.coefficients.ndim
	except KeyError:
		table.order = []
		for i in range(0,table.coefficients.ndim):
			table.order.append(data.header['ORDER%d' % i])
		order = table.order

	try:
		extents = file['EXTENTS'].data
	except KeyError:
		extents = []
	
	if len(extents) != 2*table.coefficients.ndim:
		extents = []
		for i in range(0,table.coefficients.ndim):
			extents += [table.knots[i][order[i]], table.knots[i][-order[i]-1]]

	table.extents = list(zip(extents[:-1:2], extents[1::2]))

	base_keys = set(('SIMPLE', 'BITPIX', 'EXTEND', 'TYPE'))
	for k in data.header.keys():
		if k in base_keys or any((k.startswith(prefix) for prefix in ('NAXIS', 'ORDER', 'PERIOD'))):
			continue
		setattr(table, k.lower(), data.header[k])

	file.close()
	return table

