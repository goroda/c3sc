#ifndef NODEUTIL_H
#define NODEUTIL_H

#include <stdlib.h>

int transition_assemble(size_t dx, size_t du, size_t dw, double h, double * hvec,
                        const double * drift, const double * grad_drift,
                        const double * ddiff, const double * grad_ddiff,
                        double * prob, double * grad_prob,
                        double * dt,   double * grad_dt, double * space);

int convert_fiber_to_ind(size_t d, size_t N, const double ** x, 
                         const size_t * Ngrid, double ** xgrid,
                         size_t * fixed_ind, size_t * dim_vary);
#endif
