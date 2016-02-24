#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "c3sc_elements.h"

struct State * 
euler_step(struct State *, struct Control *, double,
           struct Drift *, double *);

struct State * 
euler_maruyama_step(struct State *, double *,
                   struct Control *, double,
                   struct Dyn *,
                    double *, double *);
#endif
