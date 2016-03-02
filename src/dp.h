#ifndef DP_H
#define DP_H

#include "c3sc_elements.h"



struct DPih;

void dpih_init_ref(struct DPih *, struct Boundary *,
                   struct TensorMM *, struct Cost *,
                   struct Policy *, double,
                   int (*)(double,double*,double*,double*),
                   int (*)(double,double*,double*));

void dpih_add_transform_ref(struct DPih *,
                            struct LinTransform *, 
                            double *);

double dpih_pi_iter_approx(struct DPih *,int);

#endif
