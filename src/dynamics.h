#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "c3sc_elements.h"

void drift_init(struct Drift *,size_t,size_t,
                double*,double*,double*,double*);
int drift_eval(struct Drift *,double,double *,
               double *,double*);

////////////////////////////////////

void diff_init(struct Diff *,size_t,size_t,size_t,
               double*,double*,double*,double*);
int diff_eval(struct Diff *,double,double *,double *,double*);
size_t diff_getdw(struct Diff *);


////////////////////////////////////////////////


int dyn_eval(struct Dyn *,double,double *,double *,double *,
             double *);
size_t dyn_getdw(struct Dyn *);
    
int std_dyn_eval(struct StdDyn *, double , double *, double *,
                 double *, double *);
#endif
