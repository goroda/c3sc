#ifndef C3SC_UTIL_H
#define C3SC_UTIL_H

int c3sc_check_bounds(size_t,double*,double*,double*);
size_t c3sc_sample_discrete_rv(size_t,double *,double);
double * c3sc_combine_and_sort(size_t,double *,size_t,double *,size_t *);

#endif
