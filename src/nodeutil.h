#ifndef NODEUTIL_H
#define NODEUTIL_H

#include <stdlib.h>
#include "boundary.h"
#include "valuefunc.h"
int transition_assemble(size_t dx, size_t du, size_t dw, double h, double * hvec,
                        const double * drift, const double * grad_drift,
                        const double * ddiff, const double * grad_ddiff,
                        double * prob, double * grad_prob,
                        double * dt,   double * grad_dt, double * space);

int process_fibers(size_t d, size_t nvals, const double * x,
                   int * boundv, int * neighbors, 
                   struct Boundary * bound);

int process_fibers_neighbor(size_t d, const size_t * fixed_ind, size_t dim_vary, 
                            const double * x, int * absorbed, size_t * neighbors_vary,
                            size_t * neighbors_fixed, const size_t * ngrid, 
                            const struct Boundary * bound); 

int convert_fiber_to_ind(size_t d, size_t N, const double * x, 
                         const size_t * Ngrid, double ** xgrid,
                         size_t * fixed_ind, size_t * dim_vary);

int mca_get_neighbor_costs(size_t d,size_t N,const double * x,struct Boundary * bound,
                           struct ValueF * vf, const size_t * ngrid, double ** xgrid,
                           size_t *, size_t *,
                           int * absorbed, double * out);
#endif
