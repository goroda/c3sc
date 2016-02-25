#ifndef TENS_MARKOV_H
#define TENS_MARKOV_H

#include "c3sc_elements.h"

struct Boundary
{
    int (*bcheck)(double *, void *);
    void * args;
};

// Tensor markov model;
struct TensorMM
{
    size_t d;
    double h;
    struct Dyn * dyn;
    struct Boundary * bound;
    double * space;
};

void tensor_mm_init_ref(struct TensorMM *, size_t, 
                        double h, struct Dyn *,
                        struct Boundary *);

int tensor_mm_tprob(struct TensorMM *, double,
                    double *, double *,
                    double *, double *);


#endif
