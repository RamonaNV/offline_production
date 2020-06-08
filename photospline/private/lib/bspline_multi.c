/*
 * bspline_multi.c: Provides efficient routines using vector intrinsics to
 *    simultaneously compute the value and gradient of an N-dimensional
 *    B-spline surface.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if defined(__i386__) || defined (__x86_64__)
#ifdef __GLIBC__
#include <alloca.h>
#endif
#include <xmmintrin.h>
#elif defined(__powerpc__)
#include <altivec.h>
#endif

#include "photospline/bspline.h"

#define MAXDIM	    8
#define VECTOR_SIZE 4

#if __GNUC__ == 3
#if VECTOR_SIZE != 4
    #error On GCC 3, VECTOR_SIZE must be 4!
#endif
typedef float v4sf __attribute__(( mode(V4SF) ));
#else
typedef float v4sf __attribute__((vector_size(VECTOR_SIZE*sizeof(float))));
#endif

#define NVECS MAXDIM/VECTOR_SIZE

#if defined(__i386__) || defined (__x86_64__)
#define v4sf_init(a, b) a = _mm_set1_ps(b)
#elif defined(__powerpc__)
#ifdef vec_splats
#define v4sf_init(a, b) a = vec_splats(b)
#else
#define v4sf_init(a, b) { float b_tmp __aligned(16) = b; \
    a = vec_splat(*((v4sf *)(&b_tmp)), 0); }
#endif
#else
#define v4sf_init(a, b) { \
	((float *)(&a))[0] = b; \
	((float *)(&a))[1] = b; \
	((float *)(&a))[2] = b; \
	((float *)(&a))[3] = b; \
}
#endif

static int
maxorder(int *order, int ndim)
{
	int i, max = 0;
	
	for (i = 0; i < ndim; i++)
		if (order[i] > max)
			max = order[i];
	
	return (max);
}

static void 
ndsplineeval_multibasis_core(const struct splinetable *table, const int *centers,
    const v4sf **restrict localbasis[table->ndim], v4sf *restrict result)
{
#if (defined(__i386__) || defined (__x86_64__)) && defined(__ELF__)
	/*
	 * Workaround GCC ABI-compliance issue with SSE on x86 by
	 * forcibly realigning the stack to a 16-byte boundary.
	 */
	volatile register unsigned long sp __asm("esp");
	__asm("" : "=r"(sp));
	if (__builtin_expect(sp & 15UL, 0))
		(void)alloca(16 - (sp & 15UL));
#endif
	int i, j, k, n, tablepos;
	v4sf basis_tree[table->ndim+1][NVECS];
	int nchunks;
	int decomposedposition[table->ndim];

	tablepos = 0;
	for (n = 0; n < table->ndim; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - table->order[n])*table->strides[n];
	}
	
	for (k = 0; k < NVECS; k++) {
		v4sf_init(basis_tree[0][k], 1);
		for (n = 0; n < table->ndim; n++)
			basis_tree[n+1][k] = basis_tree[n][k]*localbasis[n][0][k];
	}
	
	nchunks = 1;
	for (n = 0; n < table->ndim - 1; n++)
		nchunks *= (table->order[n] + 1);

	n = 0;
	while (1) {
		for (i = 0; __builtin_expect(i < table->order[table->ndim-1] +
		    1, 1); i++) {
			v4sf weights;
			v4sf_init(weights, table->coefficients[tablepos + i]);
			for (k = 0; k < NVECS; k++)
				result[k] += basis_tree[table->ndim-1][k]*
				    localbasis[table->ndim-1][i][k]*weights;
		}

		if (__builtin_expect(++n == nchunks, 0))
			break;

		tablepos += table->strides[table->ndim-2];
		decomposedposition[table->ndim-2]++;

		/* Carry to higher dimensions */
		for (i = table->ndim-2;
		    decomposedposition[i] > table->order[i]; i--) {
			decomposedposition[i-1]++;
			tablepos += (table->strides[i-1]
			    - decomposedposition[i]*table->strides[i]);
			decomposedposition[i] = 0;
		}
		for (j = i; __builtin_expect(j < table->ndim-1, 1); j++)
			for (k = 0; k < NVECS; k++)
				basis_tree[j+1][k] = basis_tree[j][k]*
				    localbasis[j][decomposedposition[j]][k];
	}
}

void
bspline_nonzero(const double *knots, const unsigned nknots,
    const double x, int left, const int n,
    float *restrict values, float *restrict derivs)
{
	int i, j;
	double temp, a;
	double delta_r[n+1], delta_l[n+1];
	
	/* Special case for constant splines */
	if (n == 0) {
		values[0] = 1;
		derivs[0] = 0;
		return;
	}
	
	/*
	 * Handle the (rare) cases where x is outside the full
	 * support of the spline surface.
	 */
	assert(left >= n && left <= nknots-n-2);
	if (left == n)
		while (left >= 0 && x < knots[left])
			left--;
	else if (left == nknots-n-2)
		while (left < nknots-1 && x > knots[left+1])
			left++;
	
	/* Get the non-zero n-1th order B-splines at x */
	bsplvb(knots, x, left, 0, n, values, delta_r, delta_l);
	
	/* 
	 * Now, form the derivatives of the nth order B-splines from
	 * linear combinations of the lower-order splines.
	 *
	 * NB: bspline_deriv_nonzero() uses double-precision
	 *     temporaries, so we do the same here to ensure that
	 *     the results are identical.
	 */
	
	/* 
	 * On the last supported segment of the ith nth order spline,
	 * only the i+1th n-1th order spline is nonzero.
	 */
	temp = values[0];
	derivs[0] =  - n*temp / ((knots[left+1] - knots[left+1-n]));
	/* On the middle segments, both the ith and i+1th splines contribute. */
	for (i = 1; i < n; i++) {
		a = n*temp/((knots[left+i] - knots[left+i-n]));
		temp = values[i];
		derivs[i] = a - n*temp/(knots[left+i+1] - knots[left+i+1-n]);
	}
	/*
	 * On the first supported segment of the i+nth nth order spline,
	 * only the ith n-1th order spline is nonzero.
	 */
	derivs[n] = n*temp/((knots[left+n] - knots[left]));
	
	/* Now, continue to the non-zero nth order B-splines at x */
	bsplvb(knots, x, left, n-1, n+1, values, delta_r, delta_l);
	
	/* Rearrange for partially-supported points. */
	if ((i = n-left) > 0) {
		for (j = 0; j < left+1; j++) {
			values[j] = values[j+i]; /* Move valid splines over. */
			derivs[j] = derivs[j+i];
		}
		for ( ; j < n+1; j++)
			values[j] = derivs[j] = 0.0;
	} else if ((i = left+n+2-nknots) > 0) {
		for (j = n; j > i-1; j--) {
			values[j] = values[j-i];
			derivs[j] = derivs[j-i];
		}
		for ( ; j >= 0; j--)
			values[j] = derivs[j] = 0.0;
	}
}

/* Evaluate the spline surface and all its derivatives at x */

void
ndsplineeval_gradient(const struct splinetable *table, const double *x,
    const int *centers, double evaluates[table->ndim + 1])
{
	int n, i, j; /* , offset; */
	int maxdegree = maxorder(table->order, table->ndim) + 1;
	int nbases = table->ndim + 1;
	v4sf acc[NVECS];
	float valbasis[maxdegree];
	float gradbasis[maxdegree];
	assert(table->ndim>0);
	v4sf localbasis[table->ndim][maxdegree][NVECS];
	float *acc_ptr;
	const v4sf *localbasis_rowptr[table->ndim][maxdegree];
	const v4sf **localbasis_ptr[table->ndim];

	assert(table->ndim > 0);

	if (table->ndim+1 > MAXDIM) {
		fprintf(stderr, "Error: ndsplineeval_gradient() can only "
		    "process up to %d-dimensional tables. Adjust MAXDIM in "
		    "bspline_multi.c to change this.\n", MAXDIM-1);
		exit(1);
	}

		
	for (n = 0; n < table->ndim; n++) {

		/* 
		 * Compute the values and derivatives of the table->order[n]+1 non-zero
		 * splines at x[n], filling them into valbasis and gradbasis.
		 */
		bspline_nonzero(table->knots[n], table->nknots[n],
		    x[n], centers[n], table->order[n], valbasis, gradbasis);

		assert(table->order[n]>0);		
		for (i = 0; i <= table->order[n]; i++) {
			
			((float*)(localbasis[n][i]))[0] = valbasis[i];
			
			for (j = 1; j < table->ndim+1; j++) {
				if (j == 1+n)
					((float*)(localbasis[n][i]))[j] =
					    gradbasis[i];
				else
					((float*)(localbasis[n][i]))[j] =
					    valbasis[i];
			}
			
			localbasis_rowptr[n][i] = localbasis[n][i];
		}
		
		localbasis_ptr[n] = localbasis_rowptr[n];
	}

	acc_ptr = (float*)acc;

	for (i = 0; i < nbases; i++)
		acc_ptr[i] = 0;

	ndsplineeval_multibasis_core(table, centers, localbasis_ptr, acc);

	for (i = 0; i < nbases; i++)
		evaluates[i] = acc_ptr[i];
}
