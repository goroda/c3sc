#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "dynamics.h"

int euler_maruyama_step(double *, double *, double,
                        double, double *, double *,
                        double *,double *,
                        struct Dyn *);


struct Simulator
{
    size_t ntime;
    size_t dx;
    size_t du;
    size_t dw;
    double * t; // time
    double * x; // state
    double * u; // control
    double * w; // noise
    struct Dyn * f;

    double * drift;
    double * diff;
};

#endif
