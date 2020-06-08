#ifndef _BSPLINE_H
#define _BSPLINE_H

#include "splinetable.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Compute the value of the ith nth-order basis spline of a set
 * defined by knots at the point x.
 */

double bspline(const double *knots, double x, int i, int n);
double bspline_deriv(const double *knots, double x, int i, int n, unsigned order);

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

void bsplvb_simple(const double *knots, const unsigned nknots,
    double x, int left, int jhigh, float *biatx);
void bsplvb(const double *knots, const double x, const int left, const int jlow,
    const int jhigh, float *biatx,
    double *delta_l, double *delta_r);
void bspline_nonzero(const double *knots, const unsigned nknots,
    const double x, int left, const int n, float *values, float *derivs);
void bspline_deriv_nonzero(const double *knots, const unsigned nknots,
    const double x, const int left, const int n, float *biatx);

/*
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 */

double splineeval(const double *knots, const double *weights, int nknots, double x,
    int order, int center);

/*
 * Spline table based hypersurface evaluation. ndsplineeval() takes a spline
 * coefficient table, a vector at which to evaluate the surface, and a vector
 * indicating the evaluation centers, as for splineeval().
 *
 * tablesearchcenters() provides a method to acquire a centers vector
 * for ndsplineeval() using a binary search. Depending on how the table
 * was produced, a more efficient method may be available.
 */

int tablesearchcenters(const struct splinetable *table, const double *x, int *centers);

double ndsplineeval(const struct splinetable *table, const double *x, 
    const int *centers, int derivatives);
double ndsplineeval_linalg(const struct splinetable *table, const double *x, 
    const int *centers, int derivatives);

// slicing function to get 1-d slice with certain normalizaton properties
void
ndsplineeval_slice_coeffs(const struct splinetable *table, const double *x, const int *centers, double *results,int slice_dimension, int derivative, int area_norm);

/*
* Evaluates the spline surface, optionally differentiated to the given order
* in any dimension.
*/

double ndsplineeval_deriv(const struct splinetable *table, const double *x, 
    const int *centers, const unsigned *derivatives);

/* Evaluate a spline surface and all its derivatives at x */

void ndsplineeval_gradient(const struct splinetable *table, const double *x,
    const int *centers, double *evaluates);

/*
 * Convolve a table with the spline defined on a set of knots along a given 
 * dimension and store the spline expansion of the convolved surface in the
 * table. This will raise the order of the splines in the given dimension by
 * (n_knots - 1), i.e. convolving with an order-0 spline (a box function, 
 * defined on two knots) will raise the order of the spline surface by 1.
 */
int splinetable_convolve(struct splinetable *table, const int dim,
    const double *knots, size_t n_knots);

#ifdef __cplusplus
}
#endif

#endif /* _BSPLINE_H */
