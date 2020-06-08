/*
 * bspline.c: Routines for calculating values of B-splines and their
 *  derivatives, as well as efficient computation of the values of
 *  N-dimensional B-spline surfaces.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "photospline/bspline.h"

/*
 * Compute the value of the ith nth-order basis spline of a set
 * defined by knots at the point x.
 *
 * This is implemented using the De Boor algorithm, as outlined on
 * Wikipedia.
 */

double
bspline(const double *knots, double x, int i, int n)
{
	double result;

	if (n == 0) {
		/*
		 * Special case the 0th order case, where B-Splines
		 * are constant functions from one knot to the next.
		 */

		if (x >= knots[i] && x < knots[i+1])
			return 1.0;
		else
			return 0.0;
	}

	result = (x - knots[i])*bspline(knots, x, i, n-1) /
	    (knots[i+n] - knots[i]);
	result += (knots[i+n+1] - x)*bspline(knots, x, i+1, n-1) /
	    (knots[i+n+1] - knots[i+1]);

	return result;
}

/*
 * A brain-dead reimplementation of de Boor's BSPLVB, which generates
 * the values of the non-zero B-splines at x from the bottom up without
 * unnecessarily recalculating terms. 
 * 
 * NB: for bsplvb_simple(), bspline_nonzero(), and bspline_deriv_nonzero(),
 * `left' must be the index of the nearest fully-supported knot
 * span for splines of order n, i.e. n <= left <= nknots-n-2. For bsplvb(),
 * `left' must be the index of the nonzero 0-th order spline, i.e.
 * knots[left] <= x < knots[left+1].
 *
 * See Chapter X in: 
 * 
 * Carl de Boor. A Practical Guide to Splines, volume 27 of Applied
 *     Mathematical Sciences. Springer-Verlag, 1978.
 */

void
bsplvb_simple(const double *knots, const unsigned nknots,
    double x, int left, int degree, float *restrict biatx)
{
	assert(degree>0);
	int i, j;
	double saved, term;
	double delta_l[degree], delta_r[degree];
	
	biatx[0] = 1.0;
	
	/*
	 * Handle the (rare) cases where x is outside the full
	 * support of the spline surface.
	 */
	if (left == degree-1)
		while (left >= 0 && x < knots[left])
			left--;
	else if (left == nknots-degree-1)
		while (left < nknots-1 && x > knots[left+1])
			left++;	
	
	/* 
	 * NB: if left < degree-1 or left > nknots-degree-1,
	 * the following loop will dereference addresses ouside
	 * of knots[0:nknots]. While terms involving invalid knot
	 * indices will be discarded, it is important that `knots'
	 * have (maxdegree-1)*sizeof(double) bytes of padding
	 * before and after its valid range to prevent segfaults
	 * (see parsefitstable()).
	 */
	for (j = 0; j < degree-1; j++) {
		delta_r[j] = knots[left+j+1] - x;
		delta_l[j] = x - knots[left-j];
		
		saved = 0.0;
		
		for (i = 0; i < j+1; i++) {
			term = biatx[i] / (delta_r[i] + delta_l[j-i]);
			biatx[i] = saved + delta_r[i]*term;
			saved = delta_l[j-i]*term;
		}
		
		biatx[j+1] = saved;
	}
	
	/* 
	 * If left < (spline order), only the first (left+1)
	 * splines are valid; the remainder are utter nonsense.
	 */
	if ((i = degree-1-left) > 0) {
		for (j = 0; j < left+1; j++)
			biatx[j] = biatx[j+i]; /* Move valid splines over. */
		for ( ; j < degree; j++)
			biatx[j] = 0.0; /* The rest are zero by construction. */
	} else if ((i = left+degree+1-nknots) > 0) {
		for (j = degree-1; j > i-1; j--)
			biatx[j] = biatx[j-i];
		for ( ; j >= 0; j--)
			biatx[j] = 0.0;
	}
}

void
bsplvb(const double *knots, const double x, const int left, const int jlow,
    const int jhigh, float *restrict biatx,
    double *restrict delta_l, double *restrict delta_r)
{
	int i, j;
	double saved, term;

	if (jlow == 0)
		biatx[0] = 1.0;
	
	for (j = jlow; j < jhigh-1; j++) {
		delta_r[j] = knots[left+j+1] - x;
		delta_l[j] = x - knots[left-j];
		
		saved = 0.0;
		
		for (i = 0; i < j+1; i++) {
			term = biatx[i] / (delta_r[i] + delta_l[j-i]);
			biatx[i] = saved + delta_r[i]*term;
			saved = delta_l[j-i]*term;
		}
		
		biatx[j+1] = saved;
	}
}


void
bspline_deriv_nonzero(const double *knots, const unsigned nknots,
    const double x, int left, const int n, float *restrict biatx)
{
	/* NB: it might be tempting to use unsigned integers *left* and *n* here,
	   but indices into the knot vector may be negative (up to -order) before
	   the first fully-supported knot. */
	assert(n>0);
	int i, j;
	double temp, a;
	double delta_l[n], delta_r[n];
	
	/* Special case for constant splines */
	if (n == 0)
		return;
	
	/*
	 * Handle the (rare) cases where x is outside the full
	 * support of the spline surface.
	 */
	if (left == n)
		while (left >= 0 && x < knots[left])
			left--;
	else if (left == nknots-n-2)
		while (left < nknots-1 && x > knots[left+1])
			left++;
	
	/* Get the non-zero n-1th order B-splines at x */
	bsplvb(knots, x, left, 0 /* jlow */, n /* jhigh */,
	    biatx, delta_l, delta_r);
	
	/* 
	 * Now, form the derivatives of the nth order B-splines from
	 * linear combinations of the lower-order splines.
	 */
	
	/* 
	 * On the last supported segment of the ith nth order spline,
	 * only the i+1th n-1th order spline is nonzero.
	 */
	temp = biatx[0];
	biatx[0] =  - n*temp / ((knots[left+1] - knots[left+1-n]));
	
	/* On the middle segments, both the ith and i+1th splines contribute. */
	for (i = 1; i < n; i++) {
		a = n*temp/((knots[left+i] - knots[left+i-n]));
		temp = biatx[i];
		biatx[i] = a - n*temp/(knots[left+i+1] - knots[left+i+1-n]);
	}
	/*
	 * On the first supported segment of the i+nth nth order spline,
	 * only the ith n-1th order spline is nonzero.
	 */
	biatx[n] = n*temp/((knots[left+n] - knots[left]));

	/* Rearrange for partially-supported points. */
	if ((i = n-left) > 0) {
		for (j = 0; j < left+1; j++)
			biatx[j] = biatx[j+i]; /* Move valid splines over. */
		for ( ; j < n+1; j++)
			biatx[j] = 0.0; /* The rest are zero by construction. */
	} else if ((i = left+n+2-nknots) > 0) {
		for (j = n; j > i-1; j--)
			biatx[j] = biatx[j-i];
		for ( ; j >= 0; j--)
			biatx[j] = 0.0;
	}

}

double
bspline_deriv(const double *knots, double x, int i, int n, unsigned order)
{
	double result;

	if (n == 0) {
		/*
		 * Special case the 0th order case, where B-Splines
		 * are constant functions from one knot to the next.
		 */

		return 0.0;
	}
	
	if (order <= 1) {
		result = n * bspline(knots, x, i, n-1) / (knots[i+n] - knots[i]);
		result -= n * bspline(knots, x, i+1, n-1) / (knots[i+n+1] - knots[i+1]);
	} else {
		result = n * bspline_deriv(knots, x, i, n-1, order-1) / (knots[i+n] - knots[i]);
		result -= n * bspline_deriv(knots, x, i+1, n-1, order-1) / (knots[i+n+1] - knots[i+1]);
	}
	return result;
}

double
bspline_deriv_2(const double *knots, double x, int i, int n)
{
	double result;

	if (n <= 1) {
		/*
		 * Special case the 1st order case, where B-Splines
		 * are linear functions from one knot to the next.
		 */

		return 0.0;
	}
	
	result = bspline(knots, x, i, n-2) /
	    ((knots[i+n] - knots[i])*(knots[i+n-1] - knots[i]));
	result -= bspline(knots, x, i+1, n-2) *
	    (1./(knots[i+n] - knots[i]) + 1./(knots[i+n+1] - knots[i+1])) / 
	    (knots[i+n] - knots[i+1]);
	result += bspline(knots, x, i+2, n-2) / 
	    ((knots[i+n+1] - knots[i+1])*(knots[i+n+1] - knots[i+2]));
	
	result *= n*(n-1);
	
	return result;
}


/*
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 */

double
splineeval(const double *knots, const double *weights, int nknots, double x, int order,
    int center)
{
	double work = 0.0;
	int i;

	if (center < 0) {
		/* XXX: should be a binary search */
		for (center = 0; center+1 < nknots; center++) {
			if (x > knots[center] && x < knots[center+1])
				break;
		}
	
		if (center+1 >= nknots)
			return 0.0;
	}

	i = center - order;
	if (i < 0)
		i = 0;

	while (i < nknots-order-1 && i <= center) {
		work += weights[i]*bspline(knots, x, i, order);
		i++;
	}

	return work;
}

int
tablesearchcenters(const struct splinetable *table, const double *x, int *centers)
{
	int i, min, max;

	for (i = 0; i < table->ndim; i++) {
		
		/* Ensure we are actually inside the table. */
		if (x[i] <= table->knots[i][0] ||
		    x[i] > table->knots[i][table->nknots[i]-1])
			return (-1);
		
		/*
		 * If we're only a few knots in, take the center to be
		 * the nearest fully-supported knot.
		 */
		if (x[i] < table->knots[i][table->order[i]]) {
			centers[i] = table->order[i];
			continue;
		} else if (x[i] >= table->knots[i][table->naxes[i]]) {
			centers[i] = table->naxes[i]-1;
			continue;
		}

		min = table->order[i];
		max = table->nknots[i]-2;
		do {
			centers[i] = (max+min)/2;

			if (x[i] < table->knots[i][centers[i]])
				max = centers[i]-1;
			else
				min = centers[i]+1;
		} while (x[i] < table->knots[i][centers[i]] ||
		    x[i] >= table->knots[i][centers[i]+1]);

		/*
		 * B-splines are defined on a half-open interval. For the
		 * last point of the interval, move center one point to the
		 * left to get the limit of the sum without evaluating
		 * absent basis functions.
		 */
		if (centers[i] == table->naxes[i])
			centers[i]--;
	}


	return (0);
}

static int
maxorder(int *order, int ndim)
{
	int i, max = 0;
	
	for (i = 0; i < ndim; i++)
		if (order[i] > max)
			max = order[i];
	
	return (max);
}
   
/*
 * The N-Dimensional tensor product basis version of splineeval.
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 *
 * x is the vector at which we will evaluate the space
 */

static double
ndsplineeval_core(const struct splinetable *table, const int *centers, int maxdegree,
    float localbasis[table->ndim][maxdegree])
{
	assert(table->ndim>0);
	int i, j, n, tablepos;
	float result;
	float basis_tree[table->ndim+1];
	int nchunks;
	int decomposedposition[table->ndim];

	tablepos = 0;
	for (n = 0; n < table->ndim; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - table->order[n])*table->strides[n];
	}

	basis_tree[0] = 1;
	for (n = 0; n < table->ndim; n++)
		basis_tree[n+1] = basis_tree[n]*localbasis[n][0];
	nchunks = 1;
	for (n = 0; n < table->ndim - 1; n++)
		nchunks *= (table->order[n] + 1);

	result = 0;
	n = 0;
	while (1) {
		for (i = 0; __builtin_expect(i < table->order[table->ndim-1] +
		    1, 1); i++) {

			result += basis_tree[table->ndim-1]*
			    localbasis[table->ndim-1][i]*
			    table->coefficients[tablepos + i];
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
			basis_tree[j+1] = basis_tree[j]*
			    localbasis[j][decomposedposition[j]];
	}

	return result;
}

/* This function returns bspline coefficients along a given dimension, fixing the values+coefficients for the
other dimensions. Used to obtain a 1-d spline representation that can be easily further convolved. */

void
ndsplineeval_slice_coeffs(const struct splinetable *table, const double *x, const int *centers, double *results,int slice_dimension, int derivative, int area_norm )
{
	assert(table->ndim>0);
	int n;
	int maxdegree = maxorder(table->order, table->ndim) + 1; 
	float localbasis[table->ndim][maxdegree];

	if(slice_dimension==table->ndim-1)
	{
		printf("ERROR!!! slice dimension cannot be ndim-1 in this implementation!");
	}
	
	
	for (n = 0; n < table->ndim; n++) {
		
			
			bsplvb_simple(table->knots[n], table->nknots[n],
			    x[n], centers[n], table->order[n] + 1,
			    localbasis[n]);
	}

	


	int i, j, tablepos, slice_dimension_stride, num_coeffs;
	//double result;
	num_coeffs=table->nknots[slice_dimension]-table->order[slice_dimension]-1;

	memset(results, 0, sizeof(double)*num_coeffs);
	// temp_result is stored to teomporarily store coefficients for derivative calculation
	double temp_result[num_coeffs];

	float basis_tree[table->ndim+1]; // the last basis_tree dimension is unused .. 
	int nchunks;
	int decomposedposition[table->ndim];

	//int derivative_correction[table->ndim];


	

	nchunks = 1;
	for (n = 0; n < table->ndim - 1; n++)
	{	
		//derivative_correction[n]=0;
		
		if(n==slice_dimension)
		{
			continue;
		}
		nchunks *= (table->order[n] + 1);
		
	}



	// initialize overall table position
	tablepos = 0;
	for (n = 0; n < table->ndim; n++) {
		decomposedposition[n]=0;
		if(n!=slice_dimension)
		{
			
			tablepos += (centers[n] - table->order[n])*table->strides[n];
		}
		else
		{
			slice_dimension_stride=table->strides[n];
		}
	}
	
	
		
	// initialize local basis
	basis_tree[0] = 1;
	for (n = 0; n < table->ndim; n++)
	{	
		if(n==slice_dimension)
		{
			basis_tree[n+1] = basis_tree[n];
			continue;

		}
		basis_tree[n+1] = basis_tree[n]*localbasis[n][0];
	}

		
	double temp_base_eval=0.0;
	n=0;
	while (1) {
		
		for (i = 0; __builtin_expect(i < table->order[table->ndim-1] +
		    1, 1); i++) {

			// in contrast to ndsplineeval_core, save a value for each coefficient, separated by the slice dimension stride
			temp_base_eval=basis_tree[table->ndim-1]*localbasis[table->ndim-1][i];
			for(int nc=0; __builtin_expect(nc<num_coeffs,1) ;nc++)
			{
				results[nc] += temp_base_eval*table->coefficients[tablepos + i + nc*slice_dimension_stride];
			}

		}

		if (__builtin_expect(++n == nchunks, 0))
			break;

		// special case when the slicing dimension is ndim-2 .. skip over this dimension immediately by pushing tablepos forward by (order+1)*strides
		if(slice_dimension==table->ndim-2)
		{
			tablepos += table->strides[table->ndim-2]*(table->order[table->ndim-2]+1);
			decomposedposition[table->ndim-2]+=(table->order[table->ndim-2]+1);
		}
		else
		{
			// otherwise just add one stride
			tablepos += table->strides[table->ndim-2];
			decomposedposition[table->ndim-2]++;
		}
		
		// Carry to higher dimensions
		for (i = table->ndim-2;
		    decomposedposition[i] > table->order[i]; i--) {
		
			
			if(i==slice_dimension+1)
			{
				//printf("i=slicedimsino +1 .. \n");
				decomposedposition[i-2]++;
				tablepos += (table->strides[i-2]
			    - decomposedposition[i]*table->strides[i]);
				decomposedposition[i] = 0;
				// add one extra -1, since we want to skip the slicing dimension
				i=i-1;
			}
			else
			{
				decomposedposition[i-1]++;
				tablepos += (table->strides[i-1]
			    - decomposedposition[i]*table->strides[i]);
				decomposedposition[i] = 0;
			}
			
		}

		// stacks the tree basis up .. never include the dimension of interest, ie.e the slice dimension
		for (j = i; __builtin_expect(j < table->ndim-1, 1); j++)
		 {
			//printf("last loop index .. %d\n", j);
			if(j==slice_dimension)
			{
				basis_tree[j+1] = basis_tree[j];
			}
			else
			{
				basis_tree[j+1] = basis_tree[j]*
			    localbasis[j][decomposedposition[j]];
			}
		}

	}

		

	for(int nc=0; __builtin_expect(nc<num_coeffs,1) ;nc++)
	{

		if(derivative>0)
		{	
			double y_diff;
			double x_diff;
			int deriv_order=table->order[slice_dimension];

			if(nc==0)
			{
				// first one is special
				y_diff=deriv_order*results[0];
			}
			else
			{	
				// first form derivatives
				// requires two coefficient results usually ... so only start once we know >=2 coefficients
				// since table->order is really the degree, the new derivative order is just the degree of the original spline
				
				y_diff=deriv_order*(results[nc]-results[nc-1]);

			}
			x_diff=table->knots[slice_dimension][deriv_order+nc]-table->knots[slice_dimension][nc];

			temp_result[nc]=y_diff/((double)x_diff);

			
			if(area_norm)
			{	
				double norm_factor= ((double)deriv_order)/(   table->knots[slice_dimension][deriv_order+nc] - table->knots[slice_dimension][nc]);
				temp_result[nc]/=norm_factor;
			}

			// check again for area normalization
		}
		else
		{
			if(area_norm)
			{	
				int real_order=table->order[slice_dimension]+1;
				double norm_factor= ((double)(real_order))/ (   table->knots[slice_dimension][real_order+nc] - table->knots[slice_dimension][nc]);
				results[nc]/=norm_factor;
			}

			//printf("SLICE: abs res: %d %.10f\n", nc, results[nc]);
		}

	}

	// calculated the derivate .. copy temp_result into result
	if(derivative>0)
	{
		memcpy(results, temp_result, sizeof(temp_result));
	}
	

}


double
ndsplineeval(const struct splinetable *table, const double *x, const int *centers,
    int derivatives)
{
	assert(table->ndim>0);
	int n;
	int maxdegree = maxorder(table->order, table->ndim) + 1; 
	float localbasis[table->ndim][maxdegree];
	
	for (n = 0; n < table->ndim; n++) {
		if (derivatives & (1 << n)) {
			bspline_deriv_nonzero(table->knots[n], 
			    table->nknots[n], x[n], centers[n],
			    table->order[n], localbasis[n]);
		} else {
			bsplvb_simple(table->knots[n], table->nknots[n],
			    x[n], centers[n], table->order[n] + 1,
			    localbasis[n]);
		}
	}

	return ndsplineeval_core(table, centers, maxdegree, localbasis);
}

double
ndsplineeval_deriv(const struct splinetable *table, const double *x,
    const int *centers, const unsigned *derivatives)
{

	assert(table->ndim>0);
	int i, n;
	int maxdegree = maxorder(table->order, table->ndim) + 1; 
	float localbasis[table->ndim][maxdegree];
	
	for (n = 0; n < table->ndim; n++) {
		if (derivatives == NULL || derivatives[n] == 0) {
			bsplvb_simple(table->knots[n], table->nknots[n],
			    x[n], centers[n], table->order[n] + 1,
			    localbasis[n]);
		} else if (derivatives[n] == 1) {
			bspline_deriv_nonzero(table->knots[n], 
			    table->nknots[n], x[n], centers[n],
			    table->order[n], localbasis[n]);
		} else {
			for (i = 0; i <= table->order[n]; i++)
				localbasis[n][i] = bspline_deriv(
				    table->knots[n], x[n],
				    centers[n] - table->order[n] + i,
				    table->order[n], derivatives[n]);
		}
	}

	return ndsplineeval_core(table, centers, maxdegree, localbasis);
}

double
ndsplineeval_linalg(const struct splinetable *table, const double *x,
    const int *centers, int derivatives)
{
	assert(table->ndim>0);
	int totalcoeff, n;
	int coeffstrides[table->ndim];
	gsl_matrix_float *basis1, *basis2, *basis_elem;

	assert(table->ndim > 0);
	coeffstrides[table->ndim - 1] = totalcoeff = 1;
        for (n = table->ndim-1; n >= 0; n--) {
                totalcoeff *= (table->order[n] + 1);
                if (n > 0)
                        coeffstrides[n-1] = totalcoeff;
        }

	float basis1_data[totalcoeff], basis2_data[totalcoeff],
	    elem_data[maxorder(table->order, table->ndim) + 1];
	gsl_matrix_float b1, b2, be;
	basis1 = &b1; basis2 = &b2; basis_elem = &be;
	basis1->data = basis1_data;
	basis2->data = basis2_data;
	basis_elem->data = elem_data;

	/*
	 * Form outer product basis1 = basis2 x basis_elem, filling basis_elem
	 * every time with the non-zero basis functions on each axis and
	 * swapping basis1 and basis2 via tmp_basis.
	 */
	basis2->size1 = 1;
	basis2->size2 = 1;
	basis2->data[0] = 1.0;

	for (n = table->ndim-1; n >= 0; n--) {
		gsl_matrix_float *tmp_basis;
		if (derivatives & (1 << n)) {
			bspline_deriv_nonzero(table->knots[n], 
			    table->nknots[n], x[n], centers[n],
			    table->order[n], basis_elem->data);
		} else {
			bsplvb_simple(table->knots[n], table->nknots[n],
			    x[n], centers[n], table->order[n] + 1,
			    basis_elem->data);
		}

		basis_elem->size1 = table->order[n] + 1;
		basis_elem->size2 = 1;

		basis1->size2 = basis2->size2;
		basis1->size1 = basis_elem->size1;
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		    basis1->size1, basis1->size2, basis_elem->size2, 1.0,
		    basis_elem->data, 1, basis2->data, basis2->size2, 0,
		    basis1->data, basis1->size2);
		basis1->size2 = basis1->size1 * basis1->size2;
		basis1->size1 = 1;

		tmp_basis = basis1;
		basis1 = basis2;
		basis2 = tmp_basis;
	}

	/* Now basis1 is free, so fill it with the spline coefficients */
	int i, tablepos;
	int decomposedposition[table->ndim];
	tablepos = 0;
	for (n = 0; n < table->ndim; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - table->order[n])*table->strides[n];
	}

	for (i = 0; i < table->order[table->ndim-1] + 1; i++)
		basis1->data[i] = table->coefficients[tablepos + i];
	for (n = 1; n < coeffstrides[0] /* number of chunks */; n++) {
		tablepos += table->strides[table->ndim-2];
		decomposedposition[table->ndim-2]++;
		/* Carry to higher dimensions */
		for (i = table->ndim-2; __builtin_expect(i > 0 && 
		    decomposedposition[i] > table->order[i], 0); i--) {
			decomposedposition[i-1]++;
			tablepos += (table->strides[i-1]
			    - decomposedposition[i]*table->strides[i]);
			decomposedposition[i] = 0;
		}

		for (i = 0; i < table->order[table->ndim-1] + 1; i++)
			basis1->data[n*(table->order[table->ndim-1] + 1) + i] =
			    table->coefficients[tablepos + i];
	}

	/* Take the dot product */
	__builtin_prefetch(basis1->data);
	__builtin_prefetch(basis2->data);
	return cblas_sdot(totalcoeff, basis1->data, 1, basis2->data, 1);
}

