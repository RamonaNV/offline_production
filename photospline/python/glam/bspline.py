import numpy

def pbspline(knots, x, i, n, period):
	if n == 0:
		diff = (x - knots[i % len(knots)]) % period
		width = (knots[(i+1) % len(knots)] - knots[i % len(knots)]) % period
		if diff < width:
			return 1
		else:
			return 0

	span = (knots[(i+n) % len(knots)] - knots[i % len(knots)]) % period
	diff = (x - knots[i % len(knots)]) % period

	a = diff*pbspline(knots, x, i, n-1, period)/span

	span = (knots[(i+n+1) % len(knots)] - knots[(i+1) % len(knots)]) % period
	diff = (knots[(i+n+1) % len(knots)] - x) % period

	b = diff*pbspline(knots, x, i+1, n-1, period)/span
	return a+b

def bspline(knots, x, i, n):
	"""
	Evaluate the ith nth-order B-spline basis function at x.
	
	:param knots: a sequence of knot positions
	:param x: the position at which to evaluate the spline
	:param i: the index of the basis function
	:param n: the order of the spline
	:returns: float - the value of the spline at x
	"""
	if n == 0:
		if (x >= knots[i] and x < knots[i+1]):
			return 1.
		else:
			return 0.

	a = (x - knots[i])*bspline(knots, x, i, n-1)/(knots[i+n] - knots[i])
	b = (knots[i+n+1] - x)*bspline(knots, x, i+1, n-1)/(knots[i+n+1] - knots[i+1])
	return a+b
	
def bspline_nonzero(knots, x, n):
	"""
	Evaluate the possibly non-zero B-splines at x.
	
	:param knots: a sequence of knot positions
	:param x: the position at which to evaluate the splines
	:param n: the order of spline
	:returns: ndarray - all non-zero B-splines evaluated at x
	"""
	left = numpy.digitize([x], knots)[0] - 1
	print(len(knots), left)
	jhigh = n + 1
	biatx = numpy.zeros(jhigh)
	delta_l = numpy.zeros(jhigh)
	delta_r = numpy.zeros(jhigh)
	knits = numpy.zeros(knots.size+n)
	knits[:knots.size] = knots
	bsplvb(knits, x, left, 0, jhigh, biatx, delta_l, delta_r)
	
	# handle the partial support cases
	if (left < n):
		i = n-left
		for j in range(n-i+1):
			biatx[j] = biatx[j+i]
		biatx[n-i+1:] = 0
	elif (left > len(knots)-n-2):
		i = left+jhigh+1-len(knots)
		for j in range(n, i-1, -1):
			biatx[j] = biatx[j-i]
		biatx[:i] = 0
		
	return biatx
	
def bspline_deriv_nonzero(knots, x, n):
	"""Calculate the value of the possibly non-zero b-spline derivatives at x"""
	left = numpy.digitize([x], knots)[0] - 1
	jhigh = n
	biatx = numpy.zeros(jhigh+1)
	# special case: everything's constant
	if n == 0:
		return biatx
	delta_l = numpy.zeros(jhigh)
	delta_r = numpy.zeros(jhigh)
	# calculate the values of the n-1th order B-splines at x
	bsplvb(knots, x, left, 0, n, biatx, delta_l, delta_r)
	# post-process to form the 1st derivatives of each nth order B-spline
	temp = biatx[0]
	# biatx[0] contains the first nonzero n-1th order B-spline
	# On the last supported interval of the ith nth order
	# spline, only the i+1th order n-1th order spline is nonzero.
	biatx[0] =  - n*temp / ((knots[left+1] - knots[left+1-n]))
	# now, both the ith and i+1th n-1th order splines contribute
	for i in range(1, n):
		a = n*temp/((knots[left+i] - knots[left+i-n]))
		b = n*biatx[i]/(knots[left+i+1] - knots[left+i+1-n])
		temp = biatx[i]
		biatx[i] = a-b
	# only the ith n-1th order spline is nonzero on the first supported interval
	biatx[n] = n*temp/((knots[left+n] - knots[left]))
	
	# handle the partial support cases
	if (left < n):
		i = n-left
		for j in range(n-i+1):
			biatx[j] = biatx[j+i]
		biatx[n-i+1:] = 0
	elif (left > len(knots)-n-2):
		i = left+jhigh+1-len(knots)
		for j in range(n, i-1, -1):
			biatx[j] = biatx[j-i]
		biatx[:i] = 0
	return biatx
	
def bsplvb(knots, x, left, jlow, jhigh, biatx, delta_l, delta_r):
	"""
	Braindead reimplementation of de Boor's BSPLVB
	(an algorithm for generating the non-zero B-splines from the bottom up
	without unnecessarily re-calculating terms)
	"""
	
	if jlow == 0:
		biatx[0] = 1.0
		# biatx[left] = 1.0
	
	for j in range(jlow, jhigh - 1):
		
		delta_r[j] = knots[left+j+1] - x
		delta_l[j] = x - knots[left-j]
		# if left-j < 0:
		# 	delta_l[j] = 0
		# else:
		# 	delta_l[j] = x - knots[left-j]
		
		
		saved = 0.0
		
		for i in range(j+1):
			term = biatx[i] / (delta_r[i] + delta_l[j-i])
			biatx[i] = saved + delta_r[i]*term
			saved = delta_l[j-i]*term
		biatx[j+1] = saved

def bspline_deriv(knots, x, i, n, order=1):
	"""
	Calculate order-th derivative of a B-spline
	
	:param knots: knot vector
	:param x: point at which to evaluate
	:param i: index of spline
	:param n: order of spline
	:param order: degree of differentiation
	"""
	if n == 0:
		return 0.
	
	if order <= 1:
		a = n*bspline(knots, x, i, n-1)/(knots[i+n] - knots[i])
		b = n*bspline(knots, x, i+1, n-1)/(knots[i+n+1] - knots[i+1])
	else:
		a = n*bspline_deriv(knots, x, i, n-1, order-1)/(knots[i+n] - knots[i])
		b = n*bspline_deriv(knots, x, i+1, n-1, order-1)/(knots[i+n+1] - knots[i+1])
	
	return a-b
	
def bspline_deriv_2(knots, x, i, n):
	if n <= 1:
		return 0.

	a = bspline(knots, x, i, n-2)/((knots[i+n] - knots[i])*(knots[i+n-1] - knots[i]))
	b = bspline(knots, x, i+1, n-2)*(1/(knots[i+n] - knots[i]) + 1/(knots[i+n+1] - knots[i+1]))/(knots[i+n] - knots[i+1])
	c = bspline(knots, x, i+2, n-2)/((knots[i+n+1] - knots[i+1])*(knots[i+n+1] - knots[i+2]))
	return (n*(n-1))*(a - b + c)
	

def splinebasis(knots,order,x1,period = 0,spline = bspline):
	splinevals = []
	if period == 0:
		nsplines = len(knots)-order-1
	else:
		nsplines = len(knots)
	i = 0
	while i < nsplines:
		if period == 0:
			splinevals.append([spline(knots,x,i,order) for x in x1])
		else:
			splinevals.append([pbspline(knots,x,i,order,period) for x in x1])
		i = i+1

	return numpy.matrix(numpy.column_stack(splinevals))

def tsplinebasis(basis):
	"""Convert a B-spline basis to a T-spline basis."""
	return basis + basis.sum(axis=1).reshape((basis.shape[0],1)) - basis.cumsum(axis=1)
