#ifndef LOGSPLINEPDF
#define LOGSPLINEPDF

#include <gsl/gsl_rng.h>

#include "splinetable.h"

#ifdef __cplusplus
extern "C" {
#endif

void logsplinepdf_n_sample(double *result, int results, int burnin,
    double *coords, int dim, struct splinetable *table, int derivatives,
    double (* proposal)(void*), double (* proposal_pdf)(double, double, void*),
    void *proposal_info, const gsl_rng *rng);

void splinepdf_n_sample(double *result, int results, int burnin,
    double *coords, int dim, struct splinetable *table, int derivatives,
    double (* proposal)(void*), double (* proposal_pdf)(double, double, void*),
    void *proposal_info, const gsl_rng *rng);

#ifdef __cplusplus
}
#endif

#endif
