#ifndef DP_H
#define DP_H

#include "c3sc_elements.h"



struct DPih
{

    struct Boundary * bound;
    struct TensorMM * mm;
    struct Cost * cost;
    struct Policy * pol;

    double beta; // discount factor
    int (*stagecost)(double,double *,double *,double *);
    int (*boundcost)(double,double *,double *);

    int trans;
    struct LinTransform * lt;
    double * space;
};


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
