""" Convolve a B-spline surface with another B-spline along a single
dimension. This can be used to approximate a Gaussian smearing. See:

[1] Kyrre Stroem. On convolutions of b-splines. Journal of Computational and
Applied Mathematics, 55(1):1 - 29, 1994. Available from:
http://www.sciencedirect.com/science/article/B6TYH-45DHSSW-D/2/2ebe52e9%
a30fa525d01d3ba8fba37209, doi:DOI: 10.1016/0377-0427(94)90182-1.

[2] Lyle Ramshaw. Blossoming: A connect-the-dots approach to splines.
Technical Report 19, Digitial Equipment Corporation, 1987. Available from:
http://www.hpl.hp.com/techreports/Compaq-DEC/SRC-RR-19.pdf.
"""

import numpy
import copy

def divdiff(x, y):
	if len(y) == 1:
		return y[0]
	else:
		return (divdiff(x[:-1],y[:-1]) - divdiff(x[1:],y[1:]))/(x[-1] - x[0])

def factorial(x):
	out = 1.0
	num = x
	while num > 1:
		out *= num
		num -= 1
	return out

def convolved_blossom(x, y, z, bags):
	""" The local blossom of the convolution of the splines defined on knot
vectors x and y can be evaluated at point z via iterated divided differences
(see Stroem, Equation 13 and Lemma 9).

This is analogous to Stroem Equation 13 and Lemma 9, but with the prefactor
adapted to account for the fact that one of the splines is de-Boor normalized
and the other unit normalized. """

	# NB: we're convolving a de-Boor spline with a unit-norm spline,
	# hence (x[-1] - x[0]) rather than (x[-1] - x[0])*(y[-1] - y[0])
	# (as for two de-Boor splines)	
	scale = (x[-1] - x[0])
	fun_x = numpy.zeros(x.shape)
	for i, xi in enumerate(x):
		fun_y = numpy.zeros(y.shape)
		mask = xi + y - z > 0
		det = 1
		for bag in bags:
			det *= (xi + y - bag)
		fun_y[mask] = det[mask]
		fun_x[i] = divdiff(y, fun_y)
	return divdiff(x, fun_x)*scale

def convolution_row(x, xorder, y, yorder, rho, i):
	""" The coefficient of a spline in the convolved basis is a linear
combination of the blossoms of each spline in the un-convolved basis convolved
with the kernel spline. Here we just store the raw blossoms and multiply by
the coefficients of each un-convolved spline later.

This is analogous Stroem, Proposition 10, but with the prefactor adapted to
account for the fact that one of the splines is de-Boor normalized and the
other unit normalized. """

	k = xorder + 1 # degree of de-Boor basis spline
	q = yorder + 1 # degree of kernel spline
	
	nsplines = len(x) - xorder - 1
	
	bundle = rho[i:i+k+q+1]
	bags =  bundle[1:-1]
		
	blossoms = numpy.array([convolved_blossom(x[j:j+k+1], y, bundle[0], bags) for j in range(nsplines)])
	if k % 2 != 0:
		blossoms *= -1
	
	# NB: we're convolving a de-Boor spline with a unit-norm spline,
	# hence q!(k-1)! rather than (q-1)!(k-1)! (as for two de-Boor splines)
	facties = (factorial(q)*factorial(k-1))/factorial(k+q-1)
	
	return facties*blossoms

def pseudogauss_knots(order, sigma):
	if order == 0:
		return numpy.sqrt(3)*sigma*numpy.array([-1, 1])
	elif order == 1:
		return 2*sigma*numpy.array([-1, 0, 1])
	else:
		raise ValueError("I don't know how to construct an order-%d spline with variance %f" % (order, sigma))

def slice_axis(a, callback, slices=None, axis=0, current=0):
	"""Call a function on arbitrary 1-d slices of an array"""
	if slices is None:
		slices = [slice(None)]*a.ndim
	if current == a.ndim - 1:
		if current == axis:
			callback(a, slices)
		else:
			for i in range(a.shape[current]):
				myslice = list(slices)
				myslice[current] = i
				callback(a, myslice)
	elif current == axis:
		slice_axis(a, callback, slices, axis, current+1)
	else:
		for i in range(a.shape[current]):
			myslice = list(slices)
			myslice[current] = i
			slice_axis(a, callback, myslice, axis, current+1)

def convolve(spline, dim = -1, approx_order = 0, sigma = 100, verbose = False):
	"""This, my friends, is where the magic happens."""
	
	ndim = spline.coefficients.ndim
	if dim < 0:
		dim += ndim
	
	nsplines = spline.coefficients.shape[dim]
	order = spline.order[dim]
	
	k = order + 1
	q = approx_order + 1
	convorder = k + q - 1
	
	if verbose:
		print('order %d * order %d -> order %d' % (k-1, q-1, convorder))
		print(rho.size - convorder - 1,'splines')
	
	x = spline.knots[dim]
	y = pseudogauss_knots(approx_order, sigma)

	# Construct the convolved knot field
	rho = numpy.unique(numpy.array([xi + yi for xi in x for yi in y]))
	rho.sort()
	
	nsplines_conv = rho.size - convorder - 1
	
	# Construct a new spline table to hold the convolved surface
	newspline = copy.deepcopy(spline)
	shape = list(spline.coefficients.shape)
	shape[dim] = nsplines_conv
	newspline.coefficients = numpy.zeros(tuple(shape))
	newspline.order[dim] = convorder
	newspline.knots[dim] = rho
	newspline.extents[dim] = (rho[convorder], rho[-convorder-1])
	
	# Build up a matrix that transforms coefficients on the original knot vector
	# in the chosen dimension to coefficients on the convolved knot vector
	trafo = numpy.empty((nsplines_conv, nsplines))
	for i in range(nsplines_conv):
		trafo[i] = convolution_row(x, k-1, y, q-1, rho, i)

	# Apply the transformation to each slice of the knot field
	def convert(coefficients, slice_):
		newspline.coefficients[slice_] = numpy.dot(trafo, coefficients[slice_])	
	slice_axis(spline.coefficients, convert, axis=dim)
	
	return newspline