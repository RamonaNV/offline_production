from __future__ import print_function
import numpy
from .. import splinetable
from .bspline import *

from . import nnls as nnls_mockup
from scipy.optimize import nnls

def box(A,B):
	ea = numpy.ones((1,A.shape[1]),float)
	eb = numpy.ones((1,B.shape[1]),float)
	
	return numpy.matrix(numpy.asarray(numpy.kron(A, eb)) * \
	    numpy.asarray(numpy.kron(ea, B)))

def rho(A,B,p):
	sa = A.shape
	sb = B.shape

	newaxes = list(range(p,B.ndim)) + list(range(0,p))

	B = B.transpose(newaxes)
	nonp = numpy.prod([B.shape[i] for i in range(1,B.ndim)])
	B = numpy.reshape(B,(B.shape[0],nonp))
	
	C = numpy.asarray(A.transpose() * numpy.matrix(B))
	C = numpy.reshape(C,[A.shape[1]] + [sb[i] for i in newaxes[1:]])

	# Now invert the permutation we did to B on C
	invaxes = newaxes[:]
	for i in range(0,len(newaxes)): invaxes[newaxes[i]] = i
	C = C.transpose(invaxes)

	return C

def fit(z,w,coords,knots,order,smooth,
		periods=None,penalties={2:None},bases=None,iterations=1,monodim=None,dump=False):
	"""
	Fit gridded data to a B-spline surface.
	
	:param z: An ndarray containing the data to be fit
	:param w: An ndarray of weights for each data point. If you don't
	          want to weight the data points differently from each
	          other, use uniform weights (e.g. numpy.ones(z.shape)).
	:param coords: A list of coordinate axes giving the location of
	               the data points in each dimension. Each coordinate
	               axis must have the same length as the corresponding
	               dimension in *z*.
	:param knots: A sequence of knot vectors, each giving the location
	              of the B-spline knots in the cooresponding dimension.
	:param order: The order of the B-spline surface. If this is a single
	              integer the order will be the same for all dimensions;
	              you can also pass a list of integers to specify the
	              order to use in each dimension individually.
	:param smooth: The global smoothing parameter :math:`\lambda`. This
	               will be used only if no dimension-specific smoothing
	               is specified in *penalties*.
	:param periods: A list of floats, each giving the periodicity of
	                the corresponding dimension. If the periodicity is
	                nonzero, the fit will use a periodic B-spline
	                basis in that dimension where the ends are identical
	                by construction. This functionality is only
	                implemented in *glam*; *spglam* will silently ignore
	                the periodicity of the basis.
	:param penalties: A dictionary giving the penalty term to apply
	                  along each dimension. Each key-value pair
	                  specifies an order of penalty and the coefficient
	                  of that penalty order in each dimension. For
	                  example, to apply second-order penalties in the
	                  first three dimensions and a third-order penalty
	                  in the fourth (all with :math:`\lambda=1`), you
	                  would pass ``{2:[1,1,1,0], 3:[0,0,0,1]}``. If the
	                  coefficient list is ``None`` (as it is by default),
	                  the :math:`\lambda` specified by *smooth* will be
	                  used for all dimensions.
	:param monodim: If set to a non-negative integer, the spline surface
	                will be forced to be monotonic along the corresponding
	                dimension.
	                
	                .. note:: This involves solving a non-negative least \
	                   squares problem, and is thus significantly slower \
	                   than an unconstrained fit.
	:returns: SplineTable - the fit result
	"""
	ndim=z.ndim

	table = splinetable.SplineTable()
	table.knots = knots
	table.order = order

	order = numpy.asarray(order,dtype=int)
	if order.size == 1:
		order = order * numpy.ones(len(knots),dtype=int)
	
	# the user can pass an arbitrary linear combination of penalty orders
	penorder = [dict() for i in range(len(knots))]
	for o,coefficients in penalties.items():
		if int(o) <= 0:
			raise ValueError("Penalty order must by > 0 (not %s)" % o)
		if coefficients is None:
			# if no coefficient is specified, use the smoothness (old behavior)
			coefficients = smooth
		coefficients = numpy.asarray(coefficients,dtype=float)	
		if coefficients.size == 1:
			coefficients = coefficients * numpy.ones(len(knots),dtype=float)
		for i,coeff in enumerate(coefficients):
			penorder[i][int(o)] = coeff
	
	if periods == None:
		periods = numpy.zeros(len(knots))
	table.periods = periods

	nsplines = []
	for i in range(0,len(knots)):
		if periods[i] == 0:
			nsplines.append(len(knots[i])-order[i]-1)
		else:
			nsplines.append(len(knots[i]))		

	print("Calculating spline basis...")

	if bases == None:
		Basis = [splinebasis(knots[i],order[i],coords[i],periods[i]) for i in range(0,ndim)]
	else:
		Basis = [numpy.matrix(i) for i in bases]
	
	if monodim is not None:
		# multiplying by a lower-triangular matrix sums b-spline coefficients
		# to yield t-spline (cumulative) coefficients
		L = numpy.tril(numpy.ones((nsplines[monodim],nsplines[monodim])))
		# the b-spline coefficients are the differences between the t-spline coefficients
		Linv = numpy.linalg.inv(L)
		Basis[monodim] = numpy.dot(Basis[monodim],L)

	print("Calculating penalty matrix...")

	def calcP(nsplines, knots, dim, order, porders):
		nspl = nsplines[dim]
		knots = knots[dim]

		def divdiff(knots, m, i):
			# Calculate divided difference coefficients
			# in order to estimate derivatives.

			if m == 0:
				return numpy.asarray([1.])

			num = numpy.append(0,divdiff(knots,m-1,i+1)) - numpy.append(divdiff(knots,m-1,i),0)
			dem = (knots[i+order+1] - knots[i+m])/(order-(m-1))
			return num/dem
		
		if dim == monodim:
			L = numpy.tril(numpy.ones((nsplines[monodim],nsplines[monodim]),dtype=numpy.double))
		
		def penalty_matrix(penorder):
			D = numpy.zeros((nspl-penorder,nspl),dtype=float)
			for i in range(0, len(D)):
				D[i][i:i+penorder+1] = divdiff(knots,penorder,i)
				
			if dim == monodim:
				D = numpy.dot(D,L)
			return numpy.asmatrix(D)

		DtD = numpy.zeros((nspl,nspl))
		for porder,coeff in porders.items():
			if coeff == 0: continue
			D1 = penalty_matrix(porder)
			DtD += coeff * D1.transpose() * D1

		def prodterm(i):
			if (i == dim):
				return DtD
			else:
				return numpy.eye(nsplines[i],dtype=float)

		a = prodterm(0)
		i = 1
		while i < ndim:
			b = prodterm(i)
			a = numpy.kron(a,b)
			i = i+1

		return a

	P = calcP(nsplines,knots,0,order[0],penorder[0])
	for i in range(1,ndim):
		P = P + calcP(nsplines,knots,i,order[i],penorder[i])
	# P = smooth*P

	sidelen = numpy.product(nsplines)
	a = numpy.reshape(numpy.zeros(sidelen,float),nsplines)

	print("Reticulating splines...")

	n = 0
	while n < iterations:
		n = n+1

		F = w
		R = w*z
		for i in range(0,ndim):
			print("\tProcessing dimension",i)
			F = rho(box(Basis[i],Basis[i]),F,i)
			R = rho(Basis[i],R,i)

		Fshape = []
		for i in range(0,ndim): Fshape.extend([nsplines[i],nsplines[i]])
		F = numpy.reshape(numpy.asarray(F), Fshape)

		# Now transpose F: first the even axes, then the odd
		Fshape = list(range(0,F.ndim,2)) + list(range(1,F.ndim,2))
		F = F.transpose(Fshape)

		F = numpy.reshape(F,(sidelen,sidelen))
		r = numpy.reshape(R,(sidelen,1))

		F = F + P
		
		print("Computing iteration %d least squares solution..." % n)
		
		if dump:
			numpy.save('AtA',F)
			numpy.save('Atb',r)

		if n > 1:
			# fit for residual on further iterations
			x = numpy.reshape(a,r.shape)
			r = numpy.asarray(r - numpy.dot(F,x))
			resid = (r**2).sum()
			print('The sum of squared residuals is %e'%resid)
		
		if monodim is not None:
			print('%d negative coefficients' % (a < 0).sum())
		
		if monodim is None:
			result = numpy.linalg.lstsq(F, r)
		else:
			# result = nnls(numpy.asarray(F),numpy.asarray(r).flatten())
			# result = (nnls_mockup.nnls(F,r),)
			# result = (nnls_mockup.nnls_normal(F,r),)
			result = (nnls_mockup.nnls_normal_block(F,r),)
			#result = (nnls_mockup.nnls_normal_block3(F,r),)
			# result = (nnls_mockup.nnls_normal_block4(F,r),)
		
		coefficients = numpy.reshape(result[0],nsplines)
		
		if n == 1:
			a = coefficients
		else:
			a = a + coefficients			
		
	if monodim is not None:
		def t_to_b(t_coeff):
			b_coeff = numpy.dot(L,t_coeff)
			return b_coeff
		a = numpy.apply_along_axis(t_to_b,monodim,a)

	table.coefficients = a
	return table

def monotonize(table,monodim=0):
	"""Use the t-spline hammer to enforce monotonicity along one axis"""
	print("Futzing with t-spline basis")
	
	nsplines = []
	for i in range(0,len(table.knots)):
		if table.periods[i] == 0:
			nsplines.append(len(table.knots[i])-table.order[i]-1)
		else:
			nsplines.append(len(table.knots[i]))
	
	# multiplying by a lower-triangular matrix sums b-spline coefficients
	# to yield t-spline (cumulative) coefficients
	L = numpy.tril(numpy.ones((nsplines[monodim],nsplines[monodim])))
	# the b-spline coefficients are the differences between the t-spline coefficients
	Linv = numpy.linalg.inv(L)
	
	def futz(b_coefficients):
		# first, convert b-spline coefficients to t-spline coefficients
		coefficients = numpy.dot(Linv,b_coefficients)
		for i in range(len(coefficients)):
			a = coefficients[i]
			if a < 0:
				print('t-spline coeff %d = %e' % (i,a))
				if (i > 0) and (coefficients[i-1]+a >= 0): 
					# we're coming out of an (over-)ring; add back'erds
					coefficients[i-1] += a
					coefficients[i] = 0
				elif i+1 < len(coefficients):
					# we're going in to an (under-)ring; add forwards
					coefficients[i+1] += a
					coefficients[i] = 0

		# now, convert back to a b-spline basis
		coefficients = numpy.dot(L,coefficients)
		return coefficients
	table.coefficients = numpy.apply_along_axis(futz,monodim,table.coefficients)
	return table

def grideval(table, coords, bases=None):
	"""
	Evaluate a spline surface on a coordinate grid
	
	:param table: The spline surface to be evaluated
	:param coords: A sequence of coordinate axes, one
	               for each dimension in the spline table.
	:returns: an array of values, one for point in the tensor
	          product of the supplied axes 
    
	>>> spline.coefficients.ndim
	2
	>>> grideval(spline, [[0, 0.5, 1.0], [0.5, 1.0]])
	array([[ 1482.110672  ,  4529.49084473],
	       [ 1517.94447213,  4130.34708567],
	       [ 1506.64602055,  3425.31209966]])
	"""
	results = table.coefficients
	order = numpy.asarray(table.order,dtype=int)
	if order.size == 1:
		order = order * numpy.ones(len(table.knots),dtype=int)

	if bases == None:
		Basis = [splinebasis(table.knots[i], order[i],coords[i],
		    table.periods[i]) for i in range(0,len(table.knots))]
	else:
		Basis = bases

	for i in range(0,results.ndim):
		results = rho(Basis[i].transpose(),results,i)

	return results

