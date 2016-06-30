#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "c3.h"
#include "util.h"
#include "nodeutil.h"
#include "dynamics.h"

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

struct MCAparam{
    
    size_t dx;
    size_t du;

    size_t * ngrid;
    double ** xgrid;
    double hmin;
    double * hvec;

    double dt;
    double * grad_dt; // (du,1)

    double * prob; // (2dx+1,1)
    double * grad_prob; //(du * (2dx+1),1)
    
    double * workspace; //(du,1);
};

struct DPparam{

    // dynamics and storage space
    struct Dyn * dyn;
    double * drift;
    double * grad_drift;
    double * diff;
    double * grad_diff;

    // cost functions
    double discount; // discount factor
    int (*stagecost)(double,const double*,const double*,double*,double*);
    double * grad_stage;
    int (*boundcost)(double,const double*,double*);
    int (*obscost)(const double*,double*);
};

double bellman_control(size_t du, const double * u, double * grad_u,
                       double time, size_t dx, size_t dw, const double * x, 
                       int absorbed, double * costs,
                       struct DPparam * dp, struct MCAparam * mca)
{
    //get drift and diffusion
    int res;
    if (grad_u != NULL){
        for (size_t ii = 0; ii < du; ii++){
            grad_u[ii] = 0.0;
        }
        res = dyn_eval(dp->dyn,time,x,u,
                       dp->drift,dp->grad_drift,dp->diff,dp->grad_diff);
    }
    else{
        res = dyn_eval(dp->dyn,time,x,u,dp->drift,NULL,dp->diff,NULL);
    }
    assert (res == 0);

    // stage cost
    double val;
    double stage_cost;
        
    if (absorbed == 0){ // not absorbed
        if (grad_u != NULL){
            res = dp->stagecost(time,x,u,&stage_cost,dp->grad_stage);
            assert (res == 0);
            res = transition_assemble(dx,du,dw,mca->hmin,mca->hvec,
                                      dp->drift,dp->grad_drift,dp->diff,dp->grad_diff,mca->prob,
                                      mca->grad_prob,&(mca->dt),mca->grad_dt,mca->workspace);
            assert (res == 0);
            val = bellmanrhs(dx,du,stage_cost,dp->grad_stage,dp->discount,
                             mca->prob,mca->grad_prob,mca->dt,mca->grad_dt,
                             costs, grad_u);
        }
        else{
            res = dp->stagecost(time,x,u,&stage_cost,NULL);
            assert (res == 0);
            res = transition_assemble(dx,du,dw,mca->hmin,mca->hvec,
                                      dp->drift,NULL,dp->diff,NULL,mca->prob,
                                      NULL,&(mca->dt),NULL,mca->workspace);
            assert (res == 0);
            val = bellmanrhs(dx,du,stage_cost,NULL,dp->discount,mca->prob,NULL,mca->dt,NULL,
                             costs,NULL);
        }
    }
    else if (absorbed == 1){ // absorbed cost
        res = dp->boundcost(time,x,&val);
        assert (res == 0);
    }
    else if (absorbed == -1){
        res = dp->obscost(x,&val);
    }
    else{
        fprintf(stderr, "Unrecognized aborbed condition %d\n",absorbed);
        exit(1);
    }
    return val;
}

/* double bellman_wrapper(size_t dx, size_t du, size_t dw, size_t N, const double * x, */
/*                        const double * u, double * grad_u, */
/*                        struct DPparam * dp, struct Boundary * bound, */
/*                        struct ValueF * vf, struct MCAparam * mca) */
/* { */
/*     double time = 0.0; */
/*     double * costs = calloc_double(N * (2*dx+1)); */
/*     int * absorbed = calloc_int(N); */
/*     int res = mca_get_neighbor_costs(dx,N,x,bound,vf,mca->ngrid,mca->xgrid,absorbed,costs); */
    
/*     assert (res == 0); */
/*     for (size_t ii = 0; ii < N; ii++){ */
        
/*         double bellman_control(size_t du, const double * u, double * grad_u, */
/*                        double time, const double * x, int absorbed, double * costs, */
/*                        struct DPparam * dp, struct MCAparam * mca) */

/*     } */

/*     return 0.0; */
/* } */
