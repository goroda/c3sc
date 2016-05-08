#ifndef C3SC_UTIL_H
#define C3SC_UTIL_H

#include "c3.h"

struct ApproxArgs;
struct ApproxArgs * approx_args_init();
void approx_args_free(struct ApproxArgs *);
void approx_args_set_cross_tol(struct ApproxArgs *, double);
double approx_args_get_cross_tol(const struct ApproxArgs *);
void approx_args_set_round_tol(struct ApproxArgs *, double);
double approx_args_get_round_tol(const struct ApproxArgs *);
void approx_args_set_kickrank(struct ApproxArgs *, size_t);
size_t approx_args_get_kickrank(const struct ApproxArgs *);
void approx_args_set_maxrank(struct ApproxArgs *, size_t);
size_t approx_args_get_maxrank(const struct ApproxArgs *);

int c3sc_check_bounds(size_t,double*,double*,const double*);
size_t c3sc_sample_discrete_rv(size_t,double *,double);
double * c3sc_combine_and_sort(size_t,double *,size_t,double *,size_t *);

#endif
