/*
 * splinepdf.c: Routines to efficiently sample from a dimension of an
 *    N-dimensional spline surface. Uses a Metropolis-Hastings Markov Chain
 *    sampler.
 */

#include <math.h>
#include <gsl/gsl_rng.h>

#include "photospline/splinetable.h"
#include "photospline/bspline.h"
#include "photospline/splinepdf.h"

void logsplinepdf_n_sample(double *result, int results, int burnin,
    double *coords, int dim, struct splinetable *table, int derivatives,
    double (* proposal)(void*), double (* proposal_pdf)(double, double, void*),
    void* proposal_info, const gsl_rng *rng)
{
	int i;
	int centers[table->ndim];
	double val, lastval;
	double lastlogpdf = -INFINITY, logpdf;
	double lastproppdf, proppdf;
	double odds;
	double mint, maxt;
#ifdef DEBUG
	int accepted=0;
#endif

	/* Find the boundaries by choosing the first and last points
	 * with full support */
	mint = table->extents[dim][0];
	maxt = table->extents[dim][1];

	do {
		/*
		 * Get a starting point from the proposal distribution,
		 * making sure the PDF at the starting point is finite
		 * to avoid numerical problems.
		 */

		coords[dim] = lastval = (*proposal)(proposal_info);

		lastproppdf = (*proposal_pdf)(lastval,lastval,proposal_info);
		if (tablesearchcenters(table, coords, centers) != 0) {
			lastval = mint-1;
			continue;
		}
		lastlogpdf = ndsplineeval(table, coords, centers, 0);
		if (derivatives)
			lastlogpdf += log(ndsplineeval(table, coords, centers,
			    derivatives));
	} while (isnan(lastlogpdf) || lastval < mint || lastval > maxt);

	for (i = -burnin; i/burnin < results; i++) {
		coords[dim] = val = (*proposal)(proposal_info);

		/*
		 * If we ended up outside the table, reject the sample
		 * and try again
		 */
		if (val > maxt || val < mint ||
		    tablesearchcenters(table, coords, centers) != 0) 
			goto reject;
			
		logpdf = ndsplineeval(table, coords, centers, 0);
		if (derivatives)
			logpdf += log(ndsplineeval(table, coords, centers,
			    derivatives));
		proppdf = (*proposal_pdf)(val,lastval,proposal_info);
		odds = exp(logpdf - lastlogpdf);
		odds *= lastproppdf/proppdf;

		if (odds > 1. || gsl_rng_uniform(rng) < odds) {
			/* Accept this value */
			lastval = val;
			lastlogpdf = logpdf;
			lastproppdf = proppdf;
		#ifdef DEBUG
			accepted++;
		#endif
		}

		/*
		 * Lastval has whatever we decided on now. If the burn-in
		 * period has elapsed, write it to output array.
		 *
		 * NB: This code is used in both the accept and reject
		 * cases.
		 */

	     reject:
		if (i >= 0 && (i % burnin) == 0)
			result[i/burnin] = lastval;
	}
#ifdef DEBUG
	printf("Efficiency: %e\n", (double)(accepted)/(double)(i));
#endif
}

void splinepdf_n_sample(double *result, int results, int burnin,
    double *coords, int dim, struct splinetable *table, int derivatives,
    double (* proposal)(void*), double (* proposal_pdf)(double, double, void*),
    void *proposal_info, const gsl_rng *rng)
{
	int i;
	int centers[table->ndim];
	double val, lastval;
	double lastpdf, pdf;
	double lastproppdf, proppdf;
	double odds;
	double mint, maxt;
#ifdef DEBUG
	int accepted=0;
#endif

	/* Find the boundaries by choosing the first and last points
	 * with full support */
	mint = table->extents[dim][0];
	maxt = table->extents[dim][1];

	/* Get a fully-supported starting point from the proposal distribution. */
	do {
		coords[dim] = lastval = (*proposal)(proposal_info);
	} while (lastval < mint || lastval > maxt ||
	    tablesearchcenters(table, coords, centers) != 0);

	lastproppdf = (*proposal_pdf)(lastval,lastval,proposal_info);
	lastpdf = ndsplineeval(table, coords, centers, derivatives);

	for (i = -burnin; i/burnin < results; i++) {
		coords[dim] = val = (*proposal)(proposal_info);

		/*
		 * If we ended up outside the table, reject the sample
		 * and try again
		 */
		if (val > maxt || val < mint ||
		    tablesearchcenters(table, coords, centers) != 0) 
			goto reject;
			
		pdf = ndsplineeval(table, coords, centers, derivatives);
		proppdf = (*proposal_pdf)(val,lastval,proposal_info);
		odds = pdf/lastpdf;
		odds *= lastproppdf/proppdf;

		if (odds > 1. || gsl_rng_uniform(rng) < odds) {
			/* Accept this value */
			lastval = val;
			lastpdf = pdf;
			lastproppdf = proppdf;
			#ifdef DEBUG
				accepted++;
			#endif
		}

		/*
		 * Lastval has whatever we decided on now. If the burn-in
		 * period has elapsed, write it to output array.
		 *
		 * NB: This code is used in both the accept and reject
		 * cases.
		 */

	     reject:
		if (i >= 0 && (i % burnin) == 0)
			result[i / burnin] = lastval;
	}
	#ifdef DEBUG
		printf("Efficiency: %e\n", (double)(accepted)/(double)(i));
	#endif
}

