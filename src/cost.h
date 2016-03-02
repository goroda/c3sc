#ifndef C3SC_COST_H
#define C3SC_COST_H

#include "c3.h"

struct Cost;

struct Cost * cost_alloc(size_t d, double *, double *);

void cost_init_discrete(struct Cost *,size_t *,double **);

void cost_approx(struct Cost *,
                 double (*)(double *, void *),
                 void *, int);

struct FunctionTrain *
cost_approx_new(struct Cost *,
                double (*)(double *, void *),
                void *, int);

void cost_update_ref(struct Cost *,struct FunctionTrain *);
void cost_free(struct Cost *);


int cost_eval_neigh(struct Cost *, double, double *,
                    size_t, double[],double []);

int cost_eval(struct Cost *, double, double *, double *);

#endif
