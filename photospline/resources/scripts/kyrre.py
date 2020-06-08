#!/usr/bin/env python
from __future__ import print_function

"""
A pure-Python implementation of Kyrre Strom's algorithm for convolutions of
  B-spline defined functions with other B-spline defined functions.
  The algorithm can be found in "On convolutions of B-splines", Journal
  of Computational and Applied Mathematics, 55(1):1-29, 1994.

When run as a script, makes a few demo plots.
"""

import numpy as n
import pylab as p
import copy

from icecube.photospline.glam import bspline, glam

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

def make_a_spline(knots, order, z):
	"""Does the divided-differences-of-truncated-power-functions representation make sense?"""
	k = order + 1
	
	vals = n.zeros(z.shape)
	for i, zi in enumerate(z):
		tpf = (zi-knots)
		mask = tpf > 0
		tpf[mask] = tpf[mask]**(k-1)
		tpf[n.logical_not(mask)] = 0
		vals[i] = k*divdiff(knots, tpf)
	return vals

def cbar_simple(x, y, z, bags):
	"""
The local blossom of the convolution of the splines defined on knots x and y
can be evaluated at point z via iterated divided differences (see Stroem,
Equation 13 and Lemma 9).

This version is adapted for de Boor normalization.
	"""

	# NB: we're convolving a de-Boor spline with a unit-norm spline,
	# hence (x[-1] - x[0]) rather than (x[-1] - x[0])*(y[-1] - y[0])
	# (as for two de-Boor splines)	
	scale = (x[-1] - x[0])
	fun_x = n.zeros(x.shape)
	for i, xi in enumerate(x):
		fun_y = n.zeros(y.shape)
		mask = xi + y - z > 0
		det = 1
		for bag in bags:
			det *= (xi + y - bag)
		fun_y[mask] = det[mask]
		fun_x[i] = divdiff(y, fun_y)
	return divdiff(x, fun_x)*scale

def convolution_row(x, xorder, y, yorder, rho, i):
	"""
The coefficient of a spline in the convolved basis is a linear combination of the blossoms
of each spline in the un-convolved basis convolved with the kernel spline. Here we just store
the raw blossoms and multiply by the coefficients of each un-convolved spline later.
	"""
	k = xorder + 1 # degree of de-Boor basis spline
	q = yorder + 1 # degree of kernel spline
	
	nsplines = len(x) - xorder - 1
	
	bundle = rho[i:i+k+q+1]
	bags =  bundle[1:-1]
	
	blossoms = n.array([cbar_simple(x[j:j+k+1], y, bundle[0], bags) for j in range(nsplines)])
	if k % 2 != 0:
		blossoms *= -1
	
	# NB: we're convolving a de-Boor spline with a unit-norm spline,
	# hence q!(k-1)! rather than (q-1)!(k-1)! (as for two de-Boor splines)
	facties = (factorial(q)*factorial(k-1))/factorial(k+q-1)
	
	return facties*blossoms

	
def convolved_coefficient(x, xorder, y, yorder, rho, i):
	k = xorder + 1
	q = yorder + 1
	
	ndeg = k + q - 1
	
	bundle = rho[i:i+k+q+1]
	# xknots = x[i:i+k+1]
	xknots = x
	yknots = y
	
	# scale = (factorial(k)*factorial(q)/factorial(k+q-1))*(bundle[-1] - bundle[0])/(k+q)
	bags =  bundle[1:-1]
	cbar = 0
	
	nsplines = len(x) - xorder - 1
	
	cbar_vec = n.array([cbar_simple(xknots[j:j+k+1], yknots, bundle[0], bags) for j in range(nsplines)])
	# print cbar_vec.shape
	# cbar = n.dot(xcoefficients.flatten(), cbar_vec)

	if k % 2 != 0:
		cbar_vec *= -1
	
	# integral of a De Boor-normalized spline
	norm = (rho[i+k+q] - rho[i])/(k+q)
	
	scale = (bundle[-1] - bundle[0])/norm
	# return scale*cbar/norm
	return scale*cbar_vec
	# return scale*cbar
	
	
def blossom(x, y, z, bags):
	pass

from scipy.stats import gamma

def pandel(t,distance=50):
	lam = 71.; # meter
	tau = 671.; # ns
	x0 = 154.; # meter
	c = 0.2998; # meter/ns
	ng = 1.34;
	
	shape = distance/lam
	scale = 1./(1./tau + c/ng/x0)
		
	pdf = gamma.pdf(t,shape,scale=scale)
	return pdf

def pseudogauss_knots(order, sigma):
	if order == 0:
		return n.sqrt(3)*sigma*n.array([-1, 1])
	elif order == 1:
		return 2*sigma*n.array([-1, 0, 1])
	else:
		raise ValueError("I don't know how to construct an order-%d spline with variance %f" % (order, sigma))

def twiddle(spline, dim = -1, approx_order = 0, sigma = 100):
	
	newspline = copy.deepcopy(spline)
	x = spline.knots[dim]
	
	y = pseudogauss_knots(approx_order, sigma)
	# y = spread*n.array([-1, 0, 1])
	
	nsplines = spline.coefficients.shape[dim]
	order = spline.order[dim]	
	
	rho = n.unique(n.array([xi + yi for xi in x for yi in y]))
	rho.sort()
	
	k = spline.order[0] + 1
	q = len(y) - 1
	
	convorder = k + q - 1
	print('order %d * order %d -> order %d' % (k-1, q-1, convorder))
	
	matrix = []
	print(rho.size - convorder - 1,'splines')
	for i in range(rho.size - convorder - 1):
		matrix.append(convolution_row(x, k-1, y, q-1, rho, i))
	# matrix = scale*norm.reshape((1,norm.size))*n.array(matrix)
	matrix = n.array(matrix)

	
	shape = list(spline.coefficients.shape)
	shape[dim] = rho.size - convorder - 1
	newspline.coefficients = n.zeros(tuple(shape))	
	
	def slice_axis(a, callback, slices=None, axis=0, current=0):
		"""Call callback with arbitrary 1-d slices of an array"""
		if slices is None:
			slices = [slice(None)]*a.ndim
		# a = n.asarray(a)
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

	if dim < 0:
		axis = spline.coefficients.ndim + dim
	else:
		axis = dim
	
	# NB: need to multiply in normalization factors for each tensor-product spline
	def convert(coefficients, slice_):
		out = n.dot(matrix, coefficients[slice_])
		# print out.shape, newspline.coefficients[slice_].shape
		newspline.coefficients[slice_] = out
	
	slice_axis(spline.coefficients, convert, axis=axis)
	
	newspline.order[dim] = convorder
	newspline.knots[dim] = rho
	
	return newspline

def test_2d():
	def funky(x, y):
		return (x-1)**2 + abs(y-1)
		
	xgrid = n.linspace(0, 2, 20)
	ygrid = n.linspace(0, 2, 20)
	X, Y = n.meshgrid(xgrid, ygrid)
	
	z = funky(X, Y)
	w = n.ones(z.shape)
	
	xknots = n.concatenate(([-2, -1, -0.5], n.linspace(0, n.sqrt(2), 10)**2, [2.1, 2.2, 2.3, 2.4]))
	yknots = n.concatenate(([-2, -1, -0.5, 0], n.logspace(-2, n.log10(2), 10), [2.1, 2.2, 2.3, 2.4]))
	
	smooth = 1e-4
	
	spline = glam.fit(z, w, [xgrid, ygrid], [xknots, yknots], [2, 2], smooth, penalties = {2:[smooth]})

	raw = glam.grideval(spline, [xgrid, ygrid])

	newspline = twiddle(spline, dim = 1, approx_order = 1, sigma = 1)
	
	smoothed = glam.grideval(newspline, [xgrid, ygrid])
	
	import mpl_toolkits.mplot3d.axes3d as p3
	fig = p.figure()
	ax = p3.Axes3D(fig)
	
	ax.plot_wireframe(X, Y, raw, color='b', label='raw')
	ax.plot_wireframe(X, Y, smoothed, color='r', label='smoothed')
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	

def test_pandel():
	spline = makeyfakey()
	
	t = n.linspace(-500,2000, 500)
	delay = glam.grideval(spline, [t])
	
	p.figure()
	p.plot(t, delay, label='Raw spline (fit to Pandel)')
	
	convo = twiddle(spline, approx_order = 1, sigma = 100)
	convo_delay = glam.grideval(convo, [t])
	
	p.plot(t, convo_delay, label='Convoluted spline')
	
	print("Norm: %e -> %e" % ((delay[1:]*n.diff(t)).sum(), (convo_delay[1:]*n.diff(t)).sum()))
	
	p.legend()
	
def makeyfakey():
	
	knots = [n.concatenate(([-5, -1, -0.5], n.linspace(0, n.sqrt(2000), 15)**2, [2100,2200,2300,2400]))]
	t = n.linspace(0, n.sqrt(7000), 106)**2
	t = 0.5*(t[1:]+t[:-1])
	centers = [t]
	z = pandel(t,75)
	w = n.ones(z.shape)
	order = [2]
	smooth = 1e-4
	penalties = {2:[smooth]}
	spline = glam.fit(z, w, centers, knots, order, smooth, penalties = penalties)
	
	return spline

def test():
	
	k = 2
	
	nspl = 7
	
	x = n.logspace(-1,n.log10(k),k+nspl)
	
	x_coeffs = n.ones(nspl)
	
	# x = n.arange(k+1, dtype=float)
	y = 0.1*n.array([-2,-1, 0, 1, 2])[1:-1]
	# y = n.arange(-q/2,q/2 + 1, dtype=float)
	
	rho = n.unique(n.array([xi + yi for xi in x for yi in y]))
	rho.sort()	
	
	# k = len(x)-1
	q = len(y)-1
	
	z = n.linspace(min(x.min(),y.min()), max(x.max(),y.max())+1, 500)
	
	convorder = k + q - 1
	print('order %d * order %d -> order %d' % (k-1, q-1, convorder))
	# convolution product is a spline of degree k+q
	nsplines = rho.size - convorder - 1
	print(nsplines,'splines')
	coeffs = []
	
	# scale = (factorial(k)*factorial(q)/factorial(k+q-1))/(k+q)
	matty = []
	for i in range(nsplines):
		matty.append(convolution_row(x, k-1, y, q-1, rho, i))
		# coeffs.append(coeff)
	matty = n.array(matty)
	# coeffs = n.array(coeffs)
	# print coeffs
	coeffs = n.dot(matty, x_coeffs)
	
	xbasis = n.asarray(bspline.splinebasis(x, k-1, z))
	ybasis = n.asarray(bspline.splinebasis(y, q-1, z))
	
	basis = n.asarray(bspline.splinebasis(rho, convorder, z))
	
	evaluate = n.zeros(z.shape)
	
	for i,zi in enumerate(z):
		tpf_x = n.zeros(x.shape)
		for j,xi in enumerate(x):
			# tpf_y = xi + y - zi
			tpf_y = zi - xi + y
			mask = tpf_y > 0
			tpf_y[mask] = tpf_y[mask]**(k+q-1)
			tpf_y[n.logical_not(mask)] = 0
			tpf_x[j] = divdiff(y, tpf_y)
		
		operated = divdiff(x, tpf_x)
		evaluate[i] = (factorial(k)*factorial(q)/factorial(k+q-1))*operated
	
	p.figure()
	p.plot(z, xbasis.sum(axis=1), label='f')
	p.plot(z, ybasis.flatten(), label='g')
	
	spliff = n.dot(coeffs, basis.transpose())
		
	p.plot(z, spliff, label='f*g, spline expansion')
	
	fint = (xbasis.sum(axis=1)[1:]*n.diff(z)).sum()
	spliffint = (spliff.flatten()[1:]*n.diff(z)).sum()
	
	print("integral %f -> %f" % (fint, spliffint))
	
	p.title('Analytic spline convolution at work')
	p.legend()	

if __name__ == "__main__":
	test()
	test_pandel()
	p.show()

