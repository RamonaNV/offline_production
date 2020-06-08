import numpy as n
import scipy.linalg

def lusolve(A,b):
	lu = scipy.linalg.lu_factor(A)
	solution = scipy.linalg.lu_solve(lu,b)
	return solution
	
def cholsolve(A,b):
	chol = scipy.linalg.cho_factor(A)
	solution = scipy.linalg.cho_solve(chol,b)
	return solution


def nnls(A,b,verbose=True):
	"""A mockup of the Lawson/Hanson active set algorithm
	
	See:
	A Comparison of Block Pivoting and Interior-Point Algorithms for Linear Least Squares Problems with Nonnegative Variables
	Author(s): Luis F. Portugal, Joaquim J. Judice, Luis N. Vicente
	Source: Mathematics of Computation, Vol. 63, No. 208 (Oct., 1994), pp. 625-643
	Published by: American Mathematical Society
	Stable URL: http://www.jstor.org/stable/2153286
	"""
	# let's have the proper shape
	A = n.asarray(A)
	b = n.asarray(b).reshape(b.size)
	# step 0
	F = [] # passive (solved for) set
	G = list(range(A.shape[0])) # active (clamped to zero) set
	x = n.zeros(A.shape[0])
	
	y = -n.dot(A.transpose(),b)
	
	if verbose:
		def log(mesg):
			print(mesg)
	else:
		def log(mesg):
			pass
	
	iterations = 0
	lstsqs = 0
	while True:
		iterations += 1
		# step 1
		if len(G) == 0:
			log("Active set empty, terminating after %d iterations (%d least squares computed)" % (iterations,lstsqs))
			break # the active set is the whole set, we're done
		r_G = y[G].argmin()
		r = G[r_G]
		# print x,y
		if y[r] >= 0:
			log("Dual vector is all positive, terminating after %d iterations (%d least squares computed)" % (iterations,lstsqs))
			break # x is the optimal solution, we're done
		log("Moving %d into active set" % r)
		F.append(r); F.sort()
		G.remove(r)
		feasible = False
		while not feasible:
			# print 'F:',F
			# print 'G:',G
			# step 2
			log("Unconstrained solve: %s, %s" % (A[:,F].shape,b.shape))
			x_F = n.linalg.lstsq(A[:,F],b)[0]
			lstsqs += 1
			if (x_F >= 0).all():
				x[F] = x_F
				feasible = True
			else:
				# if the new trial solution gained a negative element
				mask = (x_F <= 0)
				theta = x[F]/(x[F] - x_F)
				
				r_F = theta[mask].argmin()
				alpha = theta[mask][r_F]
				r = n.array(F)[mask][r_F]
				x[F] = x[F] + alpha*(x_F-x[F])
				log("Moving %d to passive set" % r)
				F.remove(r)
				G.append(r); G.sort()
		# step 3
		y[:] = 0
		y[G] = n.dot(A[:,G].transpose(),(n.dot(A[:,F],x[F])-b))
	return x

def nnls_normal(AtA,Atb,verbose=True):
	"""A mockup of the Lawson/Hanson active set algorithm for pre-formulated normal equations
	
	This version starts from the unconstrained solution (which may be moderately faster)"""
	# let's have the proper shape
	AtA = n.asarray(AtA)
	Atb = n.asarray(Atb).reshape(Atb.size)
	nvar = AtA.shape[0]
	maxiter = 3*nvar
	# step 0
	F = [] # passive (solved by unconstrained least squares) set
	G = list(range(nvar)) # active (clamped to zero) set
	x = n.zeros(nvar)

	if verbose:
		def log(mesg):
			print(mesg)
	else:
		def log(mesg):
			pass

	# variant: initialize with unconstrained solution
	x = lusolve(AtA,Atb)
	mask = x < 0
	indices = n.arange(x.size)
	F = list(indices[n.logical_not(mask)])
	G = list(indices[mask])
	x[mask] = 0
	y = n.zeros(x.size)
	y[G] = n.dot(AtA[:,F][G,:],x[F]) - Atb[...,G]
	log(F)
	log(G)
	log(x)
	log(y)

	iterations = 0
	lstsqs = 0
	
	while True:
		iterations += 1
		# step 1
		if iterations > 1:
			if len(G) == 0:
				log("Active set empty, terminating after %d iterations (%d LU solves)" % (iterations,lstsqs+1))
				break # the passive set is the whole set, we're done
			r_G = y[G].argmin()
			r = G[r_G]
			# print x,y
			if y[r] >= 0:
				log("Dual vector is all positive, terminating after %d iterations (%d LU solves)" % (iterations,lstsqs+1))
				break # x is the optimal solution, we're done
			log("Moving %d into active set" % r)
			F.append(r); F.sort()
			G.remove(r)
		feasible = False
		while not feasible:
			# print 'F:',F
			# print 'G:',G
			# step 2
			
			# x_F = n.linalg.lstsq(A[:,F],b)[0]
			# select only the bits of A^T*A that apply to coefficients F
			AtA_F = AtA[:,F][F,:]
			Atb_F = Atb[...,F]
			log("Unconstrained solve: %s, %s" % (AtA_F.shape,Atb_F.shape))
			x_F = lusolve(AtA_F,Atb_F)
			lstsqs += 1
			if (x_F >= 0).all():
				x[F] = x_F
				feasible = True
			else:
				# if the new trial solution gained a negative element,
				# find the worst offending coefficient and move it back to the passive set
				mask = (x_F <= 0)
				theta = x[F]/(x[F] - x_F)
				r_F = theta[mask].argmin()
				alpha = theta[mask][r_F]
				r = n.array(F)[mask][r_F]
				x[F] = x[F] + alpha*(x_F-x[F])
				log("Moving %d to passive set" % r)
				F.remove(r)
				G.append(r); G.sort()
		# step 3
		y[:] = 0
		y[G] = n.dot(AtA[:,F][G,:],x[F]) - Atb[...,G]
	return x

def nnls_normal_block(AtA,Atb,verbose=True):
	"""A mockup of the Portugal/Judice/Vicente block-pivoting algorithm for pre-formulated normal equations
	
	See:
	A Comparison of Block Pivoting and Interior-Point Algorithms for Linear Least Squares Problems with Nonnegative Variables
	Author(s): Luis F. Portugal, Joaquim J. Judice, Luis N. Vicente
	Source: Mathematics of Computation, Vol. 63, No. 208 (Oct., 1994), pp. 625-643
	Published by: American Mathematical Society
	Stable URL: http://www.jstor.org/stable/2153286
	"""
	# let's have the proper shape
	AtA = n.asarray(AtA)
	Atb = n.asarray(Atb).reshape(Atb.size)
	nvar = AtA.shape[0]
	maxiter = 3*nvar
	
	if verbose:
		def log(mesg):
			print(mesg)
	else:
		def log(mesg):
			pass
	
	# step 0
	F = [] # passive (solved by unconstrained least squares) set
	G = list(range(nvar)) # active (clamped to zero) set
	x = n.zeros(nvar)
	y = -Atb
	
	ninf = nvar + 1 # number of infeasible coefficients
	max_trials = 10 # number of block pivots to try before resorting to Murty's method
	
	p = max_trials
	iterations = 0
	
	while iterations < maxiter:
		iterations += 1
		if (x[F] >= 0).all() and (y[G] >= 0).all():
			log('All coefficients are positive, terminating after %d iterations' % iterations)
			break
		H1 = n.array(F)[x[F] < 0]
		H2 = n.array(G)[y[G] < 0]
		current_ninf = len(H1) + len(H2)
		if current_ninf < ninf:
			ninf = current_ninf
			p = max_trials
		elif current_ninf >= ninf:
			if p >= 1:
				p -= 1
			else: # Murty's method (pick the last infeasible coordinate)
				rmax1 = max(H1)
				rmax2 = max(H2)
				if rmax1 > rmax2:
					H1 = [rmax1]; H2 = []
				else:
					H1 = []; H2 = [rmax2]
		# shuffle infeasible coefficients between sets
		log('infeasibles: %d'%ninf)
		for r in H1:
			F.remove(r); G.append(r)
		for r in H2:
			G.remove(r); F.append(r)
		F.sort(); G.sort()
		AtA_F = AtA[:,F][F,:]
		Atb_F = Atb[...,F]
		Atb_G = Atb[...,G]
		log("Unconstrained solve for %d of %d coefficients" % (len(F),nvar))
		x[F] = cholsolve(AtA_F,Atb_F)
		x[G] = 0
		y[F] = 0
		y[G] = n.dot(AtA[:,F][G,:],x[F]) - Atb_G
	if iterations == maxiter:
		print('Hooo boy, this turned out badly')
	return x
	
def nnls_normal_block3(AtA,Atb,verbose=True):
	"""A mockup of a variant Adlers BLOCK3 algorithm for pre-formulated
	normal equations (lower bounds at zero only). This algorithm always reduces
	the magnitude of the residual, and is thus finite.

	See:
	http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.32.8173
	"""
	# let's have the proper shape
	AtA = n.asarray(AtA)
	Atb = n.asarray(Atb).reshape(Atb.size)
	nvar = AtA.shape[0]
	maxiter = 3*nvar

	if verbose:
		def log(mesg):
			print(mesg)
	else:
		def log(mesg):
			pass

	# start at the unconstrained solution, orthogonally projected into the
	# feasible space
	x = cholsolve(AtA,Atb)
	mask = x < 0
	F = list(n.nonzero(n.logical_not(mask))[0]) # passive set
	G = list(n.nonzero(mask)[0]) # active set
	H1 = [] # coefficients to be constrained
	H2 = [] # coefficients to be freed

	# keep parallel passive and active sets that reflect the current
	# state of the factorization
	set_F = []
	set_G = []
	
	x[G] = 0
	y = n.zeros(nvar)
	y[G] = n.dot(AtA[:,F][G,:],x[F]) - Atb[...,G]
	
	iterations = 0
	lstsqs = 0
	residual_calcs = 0
	# Heuristic numerical-instability-avoiding tolerance from Adlers' thesis
	KKT_TOL = n.finfo(n.double).eps * nvar * 1e3
	log("Stopping tolerance: %e" % KKT_TOL)

	def modify_factor(F,G,H1,H2):
		"""In the C implementation, we would refactorize AtA here."""
		for r in H1:
			F.remove(r)
			G.append(r)
		for r in H2:
			G.remove(r)
			F.append(r)
		if len(H1) > 0: del H1[:]
		if len(H2) > 0: del H2[:]
		F.sort()
		G.sort()

	while iterations < maxiter:
		iterations += 1
		# use set_G if the column partition is out of sync with the factorization
		# state, otherwise G
		G_ = set_G if len(set_G) > 0 else G
		if len(G_) == 0:
			log("Active set empty, terminating after %d iterations (%d factorizations, %d residual calculcations)" % (iterations,lstsqs+1,residual_calcs))
			break # the passive set is the whole set, we're done
		H2 = list(n.array(G_)[y[G_] < -KKT_TOL])
		# H1 and H2 must be disjoint, so we have some extremely silly
		# conflict resolution to do here. A coefficient could have been
		# marked for constraint at the end of the inner loop and then freed
		# for having a large Lagrange multiplier here. In this case, the modification
		# is a no-op and the index should be removed from both H1 and H2.
		if (len(H1) > 0):
			# copy list to avoid b0rking iterators
			for r in list(H2):
				if r in H1:
					H1.remove(r)
					H2.remove(r)
		if len(H2) == 0:
			log("Lagrange multipliers are all non-negative, terminating after %d iterations (%d factorizations, %d residual calculations)" % (iterations,lstsqs+1,residual_calcs))
			break # x is the optimal solution, we're done
		log("\tFreeing %d coefficients (y_min = %e)" % (len(H2), y[H2].min()))
		# get ride of the temporary partitioning
		if len(set_G) > 0: del set_G[:]
		if len(set_F) > 0: del set_F[:]
		feasible = False
		while not feasible:
			log("Unconstrained solve for %d of %d coefficients" % (len(F),nvar))
			log("\tRecomputing factorization from scratch (F[%d], G[%d], H1[%d], H2[%d)"%(len(F), len(G), len(H1), len(H2)))
			# update factorization state (clearing H1 and H2)
			print('F',F)
			print('G',G)
			print('H1',H1)
			print('H2',H2)
			modify_factor(F, G, H1, H2)
			AtA_F = AtA[:,F][F,:]
			Atb_F = Atb[...,F]
			x_F = cholsolve(AtA_F,Atb_F)
			lstsqs += 1
			infeasible = x_F < 0
			if not infeasible.any():
				# the new solution doesn't violate any constraints,
				# so F doesn't change
				x[F] = x_F
				feasible = True
				log("\tStep is entirely feasible")
			elif (x[F][infeasible] < KKT_TOL).all():
				# In this case, any appreciable movement along the descent direction will encounter
				# a boundary, so we can't improve the solution at all without modifying
				# the constraints. Constrain any coefficients at the boundary.
				H1 += list(n.array(F)[infeasible])
				x[H1] = 0
				log("\tConstraining %d/%d coefficients (descent at boundary)" % (infeasible.sum(), infeasible.size))
			else:
				# the new solution violates at least one constraint => find the longest
				# distance we can move toward the new solution without increasing the
				# value of the objective function
				
				# the vector from the current solution to the minimum of ||AtA[F,:][:,F] - b[F]||
				descent = x_F - x[F]
				# Ideally, we'd like to scale the descent vector continuously and find the largest
				# one that still reduces the magnitude of the residual, but that way lies madness.
				# Instead, choose from amongst the set of descent scales that would make one of 
				# the coefficients in F exactly zero.
				alphas = (x[F]/(x[F]-x_F))[infeasible]
				# sort the candidate descent scales in descending order,
				# starting with 1 (a jump directly to the subspace minimum)
				alphas = [1.0] + list(n.sort(alphas[(alphas < 1)&(alphas > KKT_TOL)])[::-1])
				
				def subresidual(x):
					return n.dot(x[F].transpose(), n.dot(AtA_F, x[F])) - 2*n.dot(x[F].transpose(), Atb_F)
				
				residual = subresidual(x)
				residual_calcs += 1
				gotcha = False
				for alpha in alphas:
					x_candidate = x.copy()
					x_candidate[F] += alpha*descent
					# project canidate into feasible space
					candidate_infeasibles = list(n.array(F)[x_candidate[F] < 0])
					x_candidate[x_candidate < 0] = 0
					candidate_residual = subresidual(x_candidate)
					residual_calcs += 1
					# find the largest step that reduces the residual in the overall problem
					if candidate_residual <= residual:
						# we can reduce the objective function by moving towards the subspace minimum
						# update the solution and zero newly-infeasible coefficients
						log("\talpha = %- e, d_residual = %- e, residual = %- .20e" % (alpha,candidate_residual-residual,candidate_residual))
						x[F] = x_candidate[F]
						#print "xcand: " + " ".join(["%- .1e"%f for f in x_candidate[F]])
						# XXX: the partition changes here. when to update F?
						H1 += candidate_infeasibles
						#print 'H1:',H1
						H1.sort()

						x[G] = 0
						feasible = True
						gotcha = True
						break
				if not gotcha:
					# we've encountered the worst case, in which no reduction of the objective function is possible.
					# Constrain all the infeasible coefficients and try again.
					H1 += list(n.array(F)[infeasible])
					x[H1] = 0
					log("\tConstraining %d/%d coefficients (alpha[%d] = 0)" % (infeasible.sum(), infeasible.size, len(alphas)-1))
				
		if len(H1) > 0:
			# the state of the factorization and that of the set partition
			# are temporarily out of sync. Move any negative coefficients to
			# a temporary passive set for the calculation of the Lagrange
			# multipliers.
			set_F = list(F)
			set_G = list(G)
			for r in H1:
				set_F.remove(r)
				set_G.append(r)
			set_G.sort()
			F_ = set_F
			G_ = set_G
			log("\t\tSets out of sync (%d negatives in F)!"%(len(H1)))
		else:
			F_ = F
			G_ = G
		# now that we've made x[F] feasible, update constrained part
		#print 'H1',H1
		#print 'F_',F_
		#print 'G_',G_
		y[:] = 0
		x[G_] = 0
		y[G_] = n.dot(AtA[:,F_][G_,:],x[F_]) - Atb[...,G_]
		#print "x_full: " + " ".join(["%- .1e"%f for f in x])
		#print "y_full: " + " ".join(["%- .1e"%f for f in y])
	if iterations == maxiter:
		print('Hooo boy, this turned out badly')
	return x
		   
def nnls_normal_block4(AtA,Atb,verbose=True):
	"""A mockup of a variant Adlers BLOCK4 algorithm for pre-formulated
	normal equations (lower bounds at zero only). This algorithm always reduces
	the magnitude of the residual, and is thus finite.
	
	CURRENTLY QUITE BROKEN.

	See:
	http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.32.8173
	"""
	# let's have the proper shape
	AtA = n.asarray(AtA)
	Atb = n.asarray(Atb).reshape(Atb.size)
	nvar = AtA.shape[0]
	maxiter = 3*nvar

	if verbose:
		def log(mesg):
			print(mesg)
	else:
		def log(mesg):
			pass

	# start at the unconstrained solution, orthogonally projected into the
	# feasible space
	x = cholsolve(AtA,Atb)
	mask = x < 0
	F = list(n.nonzero(n.logical_not(mask))[0]) # passive set
	G = list(n.nonzero(mask)[0]) # active set
	
	x[G] = 0
	y = n.zeros(nvar)
	y[G] = n.dot(AtA[:,F][G,:],x[F]) - Atb[...,G]
	
	iterations = 0
	lstsqs = 0
	KKT_TOL = n.finfo(n.double).eps

	def modify_factor(F,G,H1,H2):
		"""In the C implementation, we would refactorize AtA here."""
		for r in H1:
			F.remove(r)
			G.append(r)
		for r in H2:
			G.remove(r)
			F.append(r)
		F.sort()
		G.sort()

	while iterations < maxiter:
		iterations += 1
		# step 1
		H1 = n.array(F)[x[F] < 0]
		H2 = n.array(G)[y[G] < -KKT_TOL]
		if (H1.size == 0 and H2.size == 0):
			log("Optimal, terminating after %d iterations (%d solves)" % (iterations,lstsqs+1))
			break

		modify_factor(F, G, H1, H2)

		AtA_F = AtA[:,F][F,:]
		Atb_F = Atb[...,F]
		log("Unconstrained solve for %d of %d coefficients" % (len(F),nvar))
		x_F = cholsolve(AtA_F,Atb_F)
		lstsqs += 1
		infeasible = x_F < 0
		feasible = False

		if not infeasible.any():
			x[F] = x_F
			log("\tStep is entirely feasible")
		elif (x[F][infeasible] == 0).all():
			# In this case, any movement along the descent direction will encounter
			# a boundary, so we can't improve the solution at all without modifying
			# the constraints.
			# Pick an infeasible coefficient and constrain it.
			r = n.array(F)[infeasible][-1]
			# XXX HACK: force this coefficient to go negative in the next iteration
			x[r] = -KKT_TOL
			log("\tConstraining %d (descent at boundary)"%r)
		else:
			# the new solution violates at least one constraint => find the longest
			# distance we can move toward the new solution without increasing the
			# value of the objective function
			
			# the vector from the current solution to the minimum of ||AtA[F,:][:,F] - b[F]||
			descent = x_F - x[F]
			# Ideally, we'd like to scale the descent vector continuously and find the largest
			# one that still reduces the magnitude of the residual, but that way lies madness.
			# Instead, choose from amongst the set of descent scales that would make one of 
			# the coefficients in F exactly zero.
			alphas = (x[F]/(x[F]-x_F))[infeasible]
			# sort the candidate descent scales in descending order,
			# starting with 1 (a jump directly to the subspace minimum)
			alphas = [1.0] + list(n.sort(alphas[(alphas < 1)&(alphas > 1e-6)])[::-1]) + [0.0]
			
			def subresidual(x):
				return n.dot(x[F].transpose(), n.dot(AtA_F, x[F])) - 2*n.dot(x[F].transpose(), Atb_F)
			residual = subresidual(x)
			gotcha = False
			for alpha in alphas:
				x_candidate = x.copy()
				x_candidate[F] += alpha*descent
				# project canidate into feasible space
				H1 = n.array(F)[x_candidate[F] < 0]
				x_candidate[x_candidate < 0] = 0
				candidate_residual = subresidual(x_candidate)
				# find the largest step that reduces the residual in the overall problem
				if candidate_residual <= residual:
					if alpha == 0:
						# we've encountered the worst case, in which no reduction of the objective function is possible.
						# bound the first infeasible coefficient and try again.
						r = n.array(F)[infeasible][-1]
						# XXX HACK: force this coefficient to go negative in the next iteration
						x[r] = -KKT_TOL
						log("\tConstraining %d (alpha = 0)"%r)
					else:
						# we can reduce the objective function by moving towards the subspace minimum
						# update the solution and zero newly-infeasible coefficients
						log("\talpha = %- e, d_residual = %- e, residual = %- .20e" % (alpha,candidate_residual-residual,candidate_residual))
						x[F] = x_candidate[F]
					break
				
		# now that we've made x[F] feasible, update constrained part
		y[:] = 0
		x[G] = 0
		y[G] = n.dot(AtA[:,F][G,:],x[F]) - Atb[...,G]
	if iterations == maxiter:
		print('Hooo boy, this turned out badly')
	return x

def test(size=5):
	import pylab as p
	# A = n.random.uniform(size=size*size,low=-1,high=1).reshape((size,size))
	A = n.eye(size)
	x_true = n.random.uniform(size=size,low=-1,high=1)
	b = n.dot(A,x_true)
	x = nnls_normal_block4(A,b)
	p.figure()
	p.plot(x_true,drawstyle='steps-post',label='true')
	p.plot(x,drawstyle='steps-post',label='fit')
	p.legend()
	p.show()
