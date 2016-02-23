#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "c3sc_elements.h"

struct State * 
euler_step(struct State *, struct Control *, double,
           struct Drift *, double *);

int euler_maruyama_step(double *, double *, double,
                        double, double *,
                        double *,double *,
                        struct Dyn *);

#endif
