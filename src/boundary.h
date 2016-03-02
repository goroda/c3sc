#ifndef C3SC_BOUNDARY_H
#define C3SC_BOUNDARY_H

struct Boundary;
struct Boundary *
boundary_alloc(size_t,int (*)(double,double *,void *,int*), void *);
void boundary_free(struct Boundary *);
int boundary_type(struct Boundary *, double, double *);

#endif
