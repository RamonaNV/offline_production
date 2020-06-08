import numpy
import itertools as it

def meshgrid_nd(*arrays,**kwargs):
	"""
	Creates coordinate tensors from an arbitrary number of one-dimensional
	coordinate vectors. With 2 arguments, the behavior is identical to
	numpy.meshgrid.
	
	This can be useful, e.g. for evaluating a function that takes ndarrays
	on an N-dimensional without resorting to numpy.vectorize or explicit
	loops.
	
	The dimensions of the returned array are transposed: the coordinate given
	by the first argument varies most quickly, followed by the second, etc.
	This matches the behavior of numpy.meshgrid. To get an array with
	dimensions in lexographical order, pass lex_order=True: ::

	  >>> x,y,z=numpy.arange(0,3),numpy.arange(4,6),numpy.arange(7,10)
	  >>> X,Y,Z=numpy.meshgrid_nd(x,y,z,lex_order=False)
	  >>> X
	  array([[[0, 1, 2],
	          [0, 1, 2]],
	         [[0, 1, 2],
	          [0, 1, 2]],
	         [[0, 1, 2],
	          [0, 1, 2]]])
	  >>> Y
	  array([[[4, 4, 4],
	          [5, 5, 5]],
	         [[4, 4, 4],
	          [5, 5, 5]],
	         [[4, 4, 4],
	          [5, 5, 5]]])
	  >>> Z
	  array([[[7, 7, 7],
	          [7, 7, 7]],
	         [[8, 8, 8],
	          [8, 8, 8]],
	         [[9, 9, 9],
	          [9, 9, 9]]])
	  >>> X,Y,Z=numpy.meshgrid_nd(x,y,z,lex_order=True)
	  >>> X
	  array([[[0, 0, 0],
	          [0, 0, 0]],
	         [[1, 1, 1],
	          [1, 1, 1]],
	         [[2, 2, 2],
	          [2, 2, 2]]])
	  >>> Y
	  array([[[4, 4, 4],
	          [5, 5, 5]],
	         [[4, 4, 4],
	          [5, 5, 5]],
	         [[4, 4, 4],
	          [5, 5, 5]]])
	  >>> Z
	  array([[[7, 8, 9],
	          [7, 8, 9]],
	         [[7, 8, 9],
	          [7, 8, 9]],
	         [[7, 8, 9],
	          [7, 8, 9]]])
	
	"""
	asarrays = list(map(numpy.asarray,arrays))
	for ar in asarrays:
		if len(ar.shape) != 1: 
			raise ValueError("arguments must be 1-d arrays")
	dims = list(map(len,asarrays))
	out = []
	nD = len(dims)
	for i,arr,dim in zip(it.count(),asarrays,dims):
		shape = [1]*nD
		shape[nD-1-i] = dim
		x = arr.reshape(*shape)
		for j,k in it.izip(range(nD),reversed(range(nD))):
			if k==i: continue
			x = x.repeat(dims[k],axis=j)
		if kwargs.get('lex_order',False): x = x.transpose()
		out.append(x)
	return tuple(out)
	
# hook it in to numpy
numpy.meshgrid_nd = meshgrid_nd

import numpy.core.numeric as _nx
from numpy.core.numeric import asarray, zeros, newaxis, outer, \
	 concatenate, isscalar, array, asanyarray
from numpy.core.fromnumeric import product, reshape
def apply_along_axes(func1d,axis,arrs,*args):
	"""
	Apply a function to 1-D slices along the given axis.

	Execute `func1d(a, *args)` where `func1d` operates on a set of 1-D arrays and `a`
	is a 1-D slice of `arr` along `axis`.

	Parameters
	----------
	func1d : function
		This function should accept 1-D arrays. It is applied to 1-D
		slices of `arr` along the specified axis.
	axis : integer
		Axis along which `arr` is sliced.
	arrs : tuple 
		tuple of input arrays. All arrays must have the same shape
	args : any
		Additional arguments to `func1d`.

	Returns
	-------
	outarr : ndarray
		The output array. The shape of `outarr` is identical to the shape of
		`arr`, except along the `axis` dimension, where the length of `outarr`
		is equal to the size of the return value of `func1d`.  If `func1d`
		returns a scalar `outarr` will have one fewer dimensions than `arr`.

	See Also
	--------
	apply_over_axis : Apply a function over 1-D slices of a single array.
	apply_over_axes : Apply a function repeatedly over multiple axes.

	"""
	arrs = list(map(asarray,arrs))
	arr = arrs[0]
	nd = arr.ndim
	if axis < 0:
		axis += nd
	if (axis >= nd):
		raise ValueError("axis must be less than arr.ndim; axis=%d, rank=%d."
			% (axis,nd))
	ind = [0]*(nd-1)
	i = zeros(nd,'O')
	indlist = list(range(nd))
	indlist.remove(axis)
	i[axis] = slice(None,None)
	outshape = asarray(arr.shape).take(indlist)
	for arr in arrs[1:]: 
		if tuple(asarray(arr.shape).take(indlist)) != tuple(outshape):
			raise ValueError("Shape of all input arrays must match in all but the selected dimension.")
	i.put(indlist, ind)
	arglist = tuple([arr[tuple(i.tolist())] for arr in arrs]) + args
	res = func1d(*arglist)
	#  if res is a number, then we have a smaller output array
	if isscalar(res):
		outarr = zeros(outshape,asarray(res).dtype)
		outarr[tuple(ind)] = res
		Ntot = product(outshape)
		k = 1
		while k < Ntot:
			# increment the index
			ind[-1] += 1
			n = -1
			while (ind[n] >= outshape[n]) and (n > (1-nd)):
				ind[n-1] += 1
				ind[n] = 0
				n -= 1
			i.put(indlist,ind)
			arglist = tuple([arr[tuple(i.tolist())] for arr in arrs]) + args
			res = func1d(*arglist)
			outarr[tuple(ind)] = res
			k += 1
		return outarr
	else:
		Ntot = product(outshape)
		holdshape = outshape
		outshape = list(arr.shape)
		outshape[axis] = len(res)
		outarr = zeros(outshape,asarray(res).dtype)
		outarr[tuple(i.tolist())] = res
		k = 1
		while k < Ntot:
			# increment the index
			ind[-1] += 1
			n = -1
			while (ind[n] >= holdshape[n]) and (n > (1-nd)):
				ind[n-1] += 1
				ind[n] = 0
				n -= 1
			i.put(indlist, ind)
			arglist = tuple([arr[tuple(i.tolist())] for arr in arrs]) + args
			res = func1d(*arglist)
			#res = func1d(arr[tuple(i.tolist())],*args)
			outarr[tuple(i.tolist())] = res
			k += 1
		return outarr

numpy.apply_along_axes = apply_along_axes	
