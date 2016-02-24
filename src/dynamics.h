#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "c3sc_elements.h"

void drift_init(struct Drift *,size_t,size_t,
                double*,double*,double*,double*);
size_t drift_getdx(struct Drift *);
int drift_eval(struct Drift *,double,double *,
               double *,double*);

////////////////////////////////////

void diff_init(struct Diff *,size_t,size_t,size_t,
               double*,double*,double*,double*);
int diff_eval(struct Diff *,double,double *,double *,double*);
size_t diff_getdw(struct Diff *);


///////////////////////////////////////////////////////////

void dyn_init_ref(struct Dyn*,struct Drift *, struct Diff *);
void dyn_add_transform_ref(struct Dyn *, struct LinTransform *,
                           double *);
size_t dyn_getdx(struct Dyn *);
size_t dyn_getdw(struct Dyn *);
int dyn_eval(struct Dyn *,double,double *,double *,double *,
             double *);



////////////////////////////////////////////////////////////

#endif
