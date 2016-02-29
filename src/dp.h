#ifndef DP_H
#define DP_H

#include "c3sc_elements.h"

struct Boundary
{
    int (*bcheck)(double *, void *);
    void * args;
};

struct DPih
{

    struct Boundary * bound;
    struct TensorMM * mm;
    struct Cost * cost;
    struct Policy * pol;

    double beta; // discount factor
    int (*stagecost)(double,double *,double *,double *);
    int (*boundcost)(double,double *,double *);
};

double dpih_pi_iter(struct DPih *);

#endif
