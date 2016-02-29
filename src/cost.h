#ifndef C3SC_COST_H
#define C3SC_COST_H

#include "c3.h"

struct Cost
{
    int discrete;
    size_t ndisc;
    struct BoundingBox * bds;
    struct FunctionTrain * cost;
};

struct Cost * cost_alloc();
void cost_init_discrete(struct Cost *,size_t,
                        double *, double *,
                        size_t, struct FunctionTrain *);
void cost_approx(struct Cost *,
                 double (*)(double *, void *),
                 void *, int);
void cost_free(struct Cost *);


int cost_eval_neigh(struct Cost *, double, double *,
                    size_t, double[],double []);

int cost_eval(struct Cost *, double, double *, double *);

#endif
