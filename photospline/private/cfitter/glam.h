/*
 * GLAM (Generalized Linear Array Model) is an algorithm for efficiently
 *   computing N-D least squares fits on a grid. More information can be
 *   found at http://www.ma.hw.ac.uk/~iain/research/GLAM.html
 */

#ifndef GLAM_H
#define GLAM_H

#include <cholmod.h>

#include "photospline/splinetable.h"
#include "splineutil.h"

void glamfit(struct ndsparse *data, double *weights, double **coords,
    struct splinetable *out, double smooth, int *order, int *penorder,
    int monodim, int verbose, cholmod_common *c);
    
void glamfit_complex(struct ndsparse *data, double *weights, double **coords,
    struct splinetable *out, int *order, cholmod_sparse* penalty,
    int monodim, int verbose, cholmod_common *c);

cholmod_sparse* add_penalty_term(long *nsplines, double *knots, int ndim,
    int dim, int order, int porder, double scale, int mono,
    cholmod_sparse *penalty, cholmod_common *c);

#endif
