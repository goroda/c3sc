#ifndef TENS_MARKOV_H
#define TENS_MARKOV_H

#include "c3sc_elements.h"
#include "cost.h"


// Tensor markov model;
struct TensorMM
{
    size_t d;
    double h;
    struct Dyn * dyn;
//    struct Boundary * bound;
    double * space;
};

void tensor_mm_init_ref(struct TensorMM *, size_t, 
                        double, struct Dyn *,
                        double *);

int tensor_mm_dyn_eval(struct TensorMM *, double,
                       double *, double *);


int tensor_mm_tprob(struct TensorMM *, double,
                    double *, double *,
                    double *, double *);



struct State *
tensor_mm_step(struct TensorMM *,
               struct State *, double,
               struct Control *, double *,
               double *);


double tensor_mm_cost(struct TensorMM *,
                      double, double *,
                      double *,
                      struct Cost *,
                      double *,
                      double *);
#endif