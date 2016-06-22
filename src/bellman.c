#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "c3.h"
#include "util.h"

/* struct DPparam */
/* { */
/*     size_t dx; */
/*     size_t du; */

/*     struct Dyn * dyn; */
/*     struct Boundary * bound; */
/*     int (*stagecost)(double,const double*,const double*,double*,double*); */
/*     int (*boundcost)(double,const double*,double*); */
/*     int (*obscost)(const double*,double*); */

/*     double discount; */

/*     size_t * N; */
/*     double ** xgrid; */
/*     double * hvec; */
/*     double hmin; */
    
/*     // space */
/*     double * drift; */
/*     double * gdrift; */
/*     double * diff; */
/*     double * gdiff; */
/*     double * xspace; */
/* }; */

/**********************************************************//**
    Evaluate the rhs of the bellman equation

    \param[in]     dx         - dimension of state space
    \param[in]     du         - dimension of control space
    \param[in]     stage_cost - stage cost
    \param[in]     stage_grad - gradient of stage cost with respect to control
    \param[in]     discount   - discount factor
    \param[in]     prob       - (2dx+1,) transition probabilities
    \param[in]     prob_grad  - gradient of the transition probabilities
    \param[in]     dt         - time step
    \param[in]     dtgrad     - gradient of the time step
    \param[in]     cost       - (2dx+1,) cost function values of neighbors
    \param[in,out] grad       - gradient of output with respect to control

    \return value

    \note
    The probabilities are ordered as ("left", "right") for each dimension
    followed by the transition probability to one-self

**************************************************************/
double bellmanrhs(size_t dx, size_t du, double stage_cost, const double * stage_grad, 
                  double discount, const double * prob, const double * prob_grad, 
                  double dt, const double * dtgrad, const double * cost, double * grad)
{
    double ebt = exp(-discount * dt);
    double cost_to_go = cblas_ddot(2*dx+1,prob,1,cost,1);

    double out = dt * stage_cost + ebt * cost_to_go;

    if (grad != NULL){
        for (size_t jj = 0; jj < du; jj++){
            // derivative of stage cost
            grad[jj] = stage_grad[jj] * dt + dtgrad[jj]*stage_cost ;
            grad[jj] += (-discount) * dtgrad[jj] * ebt * cost_to_go;
            for (size_t ii = 0; ii < 2*dx+1; ii++){
                grad[jj] += ebt * prob_grad[jj*(2*dx+1) +ii] * cost[ii];
            }
        }
    }

    return out;
}

