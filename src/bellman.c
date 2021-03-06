// Copyright (c) 2015-2016, Massachusetts Institute of Technology
// Copyright (c) 2018, University of Michigan
//
// This file is part of the C3 for Stochastic Optimal Control (C3SC) toolbox
// Author: Alex A. Gorodetsky 
// Contact: goroda@umich.edu
// Website: https://www.alexgorodetsky.com

// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice, 
//    this list of conditions and the following disclaimer in the documentation 
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its contributors 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//Code


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* #define NDEBUG */
#include <assert.h>
#include <math.h>
#include <time.h>


#include "c3/array.h"
#include "util.h"
#include "nodeutil.h"
#include "dynamics.h"
#include "hashgrid.h"
#include "bellman.h"

#ifdef _OPENMP
#include "omp.h"
#endif

struct Memory
{
    void * shared;
    size_t private;
};


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
    assert (prob != NULL);

    double ebt = exp(-discount * dt);
    double cost_to_go = cblas_ddot(2*dx+1,prob,1,cost,1);

    double out = dt * stage_cost + ebt * cost_to_go;

    if (grad != NULL){
        for (size_t jj = 0; jj < du; jj++){
            // derivative of stage cost
            grad[jj] = stage_grad[jj] * dt + dtgrad[jj]*stage_cost ;
            grad[jj] += (-discount) * dtgrad[jj] * ebt * cost_to_go;
            for (size_t ii = 0; ii < 2*dx+1; ii++){
                grad[jj] += ebt * prob_grad[ii*du + jj] * cost[ii];
            }
        }
        /* printf("grad = ");dprint(du,grad); */
    }

    return out;
}


////////////////////////////////////////////////
/// MCA
///////////////////////////////////////////////
struct MCAparam{
    
    size_t dx;
    size_t du;

    size_t * ngrid;
    double ** xgrid;
    double hmin;
    double * hvec;


    double h2;
    double * t;
    
};

/**********************************************************//**
    Create/allocate memory for the MCA structure that stores 
    all the information
    needed for Kushner's solution method

    \param[in]     dx         - dimension of state space
    \param[in]     du         - dimension of control space

    \return allocated structure
**************************************************************/
struct MCAparam * mca_param_create(size_t dx, size_t du)
{
    struct MCAparam * mca = malloc(sizeof(struct MCAparam));
    assert (mca != NULL);
    mca->dx = dx;
    mca->du = du;

    mca->ngrid = NULL;
    mca->xgrid = NULL;
    mca->hvec = NULL;
    mca->hmin = 0;

    mca->t = NULL;
    
    return mca;
}


/**********************************************************//**
    Add references to the grid information to the MCA

    \param[in,out] mca   - MCA structure to modify    
    \param[in]     ngrid - size of discretization in each dimensiona (mca->dx,) array
    \param[in]     xgrid - nodes of discretization in each dimension
    \param[in]     hmin  - minimimum discretization level
    \param[in]     hvec  - distance between nodes in each dimension (mca->dx,)
**************************************************************/
void mca_add_grid_refs(struct MCAparam * mca, size_t * ngrid, double ** xgrid,
                       double hmin, double * hvec)
{
    assert (mca != NULL);
    mca->ngrid = ngrid;
    mca->xgrid = xgrid;
    mca->hmin = hmin;
    mca->hvec = hvec;


    mca->h2 = mca->hmin * mca->hmin;
    mca->t = calloc_double(2*mca->dx);
    for (size_t ii = 0; ii < mca->dx; ii++){
        mca->t[2*ii] = mca->h2 / mca->hvec[ii];
        mca->t[2*ii+1] = mca->t[2*ii] / mca->hvec[ii];
    }

}

/**********************************************************//**
    Destroy memory allocated to MCA structure
                                                           
    \param[in,out] mca   - MCA structure to modify    
**************************************************************/
void mca_param_destroy(struct MCAparam * mca)
{
    if (mca != NULL){
        free(mca->t);mca->t = NULL;
        free(mca); mca = NULL;
    }
}

////////////////////////////////////////////////
/// DP
///////////////////////////////////////////////
struct DPparam{

    // dynamics and storage space
    struct Drift * dyn_drift;
    struct Diff * dyn_diff;
    struct Boundary * bound;
    
    // cost functions
    double discount; // discount factor
    int (*stagecost)(double,const double*,const double*,double*,double*);
    int (*boundcost)(double,const double*,double*);
    int (*obscost)(const double*,double*);
};

struct DPparam * dp_param_create(size_t dx, size_t du, size_t dw, double discount)
{
    struct DPparam * dp = malloc(sizeof(struct DPparam));
    assert (dp != NULL);
    
    dp->dyn_drift = drift_alloc(dx,du);
    dp->dyn_diff = diff_alloc(dx,du,dw);
    dp->bound = NULL;

    dp->discount = discount;
    dp->stagecost = NULL;
    dp->boundcost = NULL;
    dp->obscost = NULL;

    return dp;
}

void dp_param_destroy(struct DPparam * dp)
{
    if (dp != NULL){
        drift_free(dp->dyn_drift); dp->dyn_drift = NULL;
        diff_free(dp->dyn_diff); dp->dyn_diff = NULL;
        free(dp); dp = NULL;
    }
}

void dp_param_add_drift(struct DPparam * dp, int (*b)(double,const double*,const double*,
                                                      double*,double*,void*),
                        void * bargs)
{
    assert (dp != NULL);
    drift_add_func(dp->dyn_drift,b,bargs);
}

void dp_param_add_diff(struct DPparam * dp, int (*s)(double,const double*,const double*,
                                                     double*,double*,void*),
                        void * sargs)
{
    assert (dp != NULL);
    diff_add_func(dp->dyn_diff,s,sargs);
}

void dp_param_add_boundary(struct DPparam * dp, struct Boundary * bound)
{
    assert (dp != NULL);
    dp->bound = bound;
}

void dp_param_add_stagecost(struct DPparam * dp, int (*stagecost)(double,const double*,
                                                                  const double*,double*,double*))
{
    assert (dp != NULL);
    dp->stagecost = stagecost;
}

void dp_param_add_boundcost(struct DPparam * dp, int (*boundcost)(double,const double*,double*))
{
    assert (dp != NULL);
    dp->boundcost = boundcost;
}

void dp_param_add_obscost(struct DPparam * dp, int (*obscost)(const double*,double*))
{
    assert (dp != NULL);
    dp->obscost = obscost;
}


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
struct ControlParams
{
    double time;
    size_t dx;
    size_t dw;

    size_t N;
    const double * x;
    
    struct DPparam * dp;
    struct MCAparam * mca;
    struct Workspace * work;
    struct c3Opt * opt;

    int res_last_grad; // result of the last gradient evaluation
    
};

struct ControlParams * control_params_create(size_t dx, size_t dw,
                                             struct DPparam * dp,
                                             struct MCAparam * mca,
                                             struct Workspace * work,
                                             struct c3Opt * opt)
{
    struct ControlParams * c = malloc(sizeof(struct ControlParams));
    assert (c != NULL);
    c->dx = dx;
    c->dw = dw;
    c->dp = dp;
    c->mca = mca;
    c->work = work;
    c->opt = opt;
    c->res_last_grad = 0;
    return c;
}

void control_params_add_time_and_states(struct ControlParams * cp, double time, size_t N, const double * x)
{
    assert (cp != NULL);
    cp->time = time;
    cp->N = N;
    cp->x = x;
}

int control_params_get_last_res(const struct ControlParams * cp)
{
    assert (cp != NULL);
    return cp->res_last_grad;
}

void control_params_destroy(struct ControlParams * cp)
{
    if (cp != NULL){
        free(cp); cp = NULL;
    }

}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**********************************************************//**
    Evaluate Bellmans equation for a certain control.

    \param[in] du     - number of control variables
    \param[in] u      - control
    \param[in] grad_u - control gradient
    \param[in] args   - pointer to problem parameters (ControlParams)

    \return value of bellmans equation
**************************************************************/
double bellman_control(size_t du, const double * u, double * grad_u, void * args)
{
    struct Memory * mem = args;
    struct ControlParams * param = mem->shared;
    size_t mem_node = mem->private;

    // unpack
    struct DPparam * dp = param->dp;
    struct MCAparam * mca = param->mca;
    struct Workspace * work = param->work;

    assert (dp != NULL);
    assert (mca != NULL);
    assert (work != NULL);

    double time = param->time;
    size_t dx = param->dx;
    size_t dw = param->dw;

    // stage cost
    double val;
    double stage_cost;

    const double * x = param->x + mem_node * dx;

    int absorbed        = *(workspace_get_absorbed        (work, mem_node));
    double * costs      = workspace_get_costs             (work, mem_node);
    double * drift      = workspace_get_drift             (work, mem_node);
    double * grad_drift = workspace_get_grad_drift        (work, mem_node);
    double * diff       = workspace_get_diff              (work, mem_node);
    double * grad_diff  = workspace_get_grad_diff         (work, mem_node);
    double * dt         = workspace_get_dt                (work, mem_node);
    double * grad_dt    = workspace_get_grad_dt           (work, mem_node);
    double * prob       = workspace_get_prob              (work, mem_node);
    double * grad_prob  = workspace_get_grad_prob         (work, mem_node);
    double * grad_stage = workspace_get_grad_stage        (work, mem_node);
    double * workspace  = workspace_get_control_size_extra(work, mem_node);
    

    if (absorbed == 0){ // not absorbed

        //get drift and diffusion
        int res;
        if (grad_u != NULL){
            for (size_t ii = 0; ii < du; ii++){
                grad_u[ii] = 0.0;
            }
            res = drift_eval(dp->dyn_drift,time,x,u,drift,grad_drift);
            assert (res == 0);
            res = diff_eval(dp->dyn_diff,time,x,u,diff,grad_diff);

        }
        else{
            /* printf("x = ");dprint(dx,x); */
            res = drift_eval(dp->dyn_drift,time,x,u,drift,NULL);
            assert (res == 0);
            res = diff_eval(dp->dyn_diff,time,x,u,diff,NULL);
        }
        assert (res == 0);
        
        if (grad_u != NULL){
            res = dp->stagecost(time,x,u,&stage_cost,grad_stage);
            assert (res == 0);
            /* res = transition_assemble(dx,du,dw,mca->hmin,mca->hvec, */
            /*                           drift,grad_drift,diff,grad_diff,prob, */
            /*                           grad_prob,dt,grad_dt,workspace); */
            res = transition_assemble(dx,du,dw,mca->h2,mca->t,
                                      drift,grad_drift,diff,grad_diff,prob,
                                      grad_prob,dt,grad_dt,workspace);

            param->res_last_grad = res;
            /* assert (res == 0); */
            val = bellmanrhs(dx,du,stage_cost,grad_stage,dp->discount,
                             prob,grad_prob,*dt,grad_dt,
                             costs, grad_u);
        }
        else{
            res = dp->stagecost(time,x,u,&stage_cost,NULL);
            assert (res == 0);
            /* res = transition_assemble(dx,du,dw,mca->hmin,mca->hvec, */
            /*                           drift,NULL,diff,NULL,prob, */
            /*                           NULL,dt,NULL,NULL); */
            res = transition_assemble(dx,du,dw,mca->h2,mca->t,
                                      drift,NULL,diff,NULL,prob,
                                      NULL,dt,NULL,NULL);
            assert (res == 0);
            val = bellmanrhs(dx,du,stage_cost,NULL,dp->discount,prob,NULL,*dt,NULL,
                             costs,NULL);
        }
    }
    else if (absorbed == 1){ // absorbed cost
        int res = dp->boundcost(time,x,&val);
        if (grad_u != NULL){
            for (size_t ii = 0; ii < du; ii++){
                grad_u[ii] = 0.0;
            }
        }
        assert (res == 0);
    }
    else if (absorbed == -1){
        int res = dp->obscost(x,&val);
        if (grad_u != NULL){
            for (size_t ii = 0; ii < du; ii++){
                grad_u[ii] = 0.0;
            }
        }
        assert (res == 0);
    }
    else{
        fprintf(stderr, "Unrecognized aborbed condition %d\n",absorbed);
        exit(1);
    }
    return val;
}


/* static void diagnose_opt_failure(struct ControlParams * param, double * u, size_t du) */
/* { */

/*     (void)(param); */
/*     fprintf(stderr,"u = "); */
/*     for (size_t ii = 0; ii < du; ii++){ */
/*         fprintf(stderr,"%3.8G ",u[ii]); */
/*     } */
/*     fprintf(stderr,"\n"); */
/* } */

/**********************************************************//**
    Find the optimal control

    \param[in]     du       - number of control variables
    \param[in,out] u        - control
    \param[in]     val      - value of bellman function at optimal control
    \param[in]     arg      - pointer to memory

    \return 0 if successful, otherwise not
**************************************************************/
int bellman_optimal(size_t du, double * u, double * val, void * arg)
{
    assert (arg != NULL);
    struct Memory * mem = arg;
    struct ControlParams * param = mem->shared;
    
    size_t mem_node = mem->private;
    int absorbed = *(workspace_get_absorbed(param->work, mem_node));
    /* printf("asborbed = %d\n",absorbed); */
    if (absorbed == 1){ // absorbed cost
        /* printf("asborbed on boundary!\n"); */
        double time = param->time;
        struct DPparam * dp = param->dp;
        const double * x = param->x + mem_node * param->dx;
        int res = dp->boundcost(time,x,val);
        for (size_t ii = 0; ii < du; ii++){
            u[ii] = 0.0;
        }
        return res;
    }
    else if (absorbed == -1){
        struct DPparam * dp = param->dp;
        const double * x = param->x + mem_node * param->dx;
        int res = dp->obscost(x,val);
        for (size_t ii = 0; ii < du; ii++){
            u[ii] = 0.0;
        }
        return res;
    }
    

    struct c3Opt * opt = c3opt_copy(param->opt);
    assert (opt != NULL);
    c3opt_add_objective(opt,&bellman_control,mem);

    if (c3opt_is_bruteforce(opt) == 1){
        c3opt_minimize(opt,u,val);
        c3opt_free(opt); opt = NULL;
        return 0;
    }

    double * lbu = c3opt_get_lb(opt);
    double * ubu = c3opt_get_ub(opt);
    double * umin = calloc_double(du);
    double * ucurr = calloc_double(du);
    double minval = 0.0;
    size_t nrand = 10; // note this is changed if du = 1 or du = 3;
    size_t npert = 2;
    double valtemp;
    int justrand = -1;
    /* int justrand = 1;     */

    // first do zero;
    /* valtemp = bellman_control(du,ucurr,NULL,arg); */
    /* c3opt_minimize(opt,ucurr,&valtemp); */
    /* minval = valtemp; */
    /* memmove(umin,ucurr,du*sizeof(double)); */
    for (size_t ii = 0; ii < du; ii++){
        u[ii] = (lbu[ii]+ubu[ii])/2;
    }
    int res = c3opt_minimize(opt,u,&valtemp);
    (void)(res);
    /* if (res < 0){ */
    /*     fprintf(stderr,"Optimization result %d\n",res); */
    /*     fprintf(stderr,"\t in bellman_optimal\n"); */
    /*     dprint(du,u); */
    /*     fprintf(stderr,"\n\n\n Rerunning with full output\n"); */

    /*     for (size_t ii = 0; ii < du; ii++){ */
    /*         u[ii] = (lbu[ii]+ubu[ii])/2; */
    /*     } */
        
    /*     c3opt_set_verbose(opt,5); */
    /*     int res2 = c3opt_minimize(opt,u,&valtemp); */
    /*     fprintf(stderr,"New optimization result %d\n",res2); */

        
    /*     diagnose_opt_failure(param,u,du); */
    /*     exit(1); */
    /* } */
    /* assert (res <= 0); */
    minval = valtemp;
    memmove(umin,u,du*sizeof(double));
    
    if (justrand == 1){
        // now compare with random samples
        /* double valtemp = 123456; */
        /* nrand = 10; */
        nrand = 1;        
        for (size_t jj = 0; jj < nrand; jj++){
            for (size_t kk = 0; kk < du; kk++){
                ucurr[kk] = randu()*(ubu[kk]-lbu[kk]) + lbu[kk];
            }
            /* valtemp = bellman_control(du,ucurr,NULL,arg); */
            c3opt_minimize(opt,ucurr,&valtemp);            
            if (valtemp < minval){
                minval = valtemp;
                memmove(umin,ucurr,du*sizeof(double));
            }
        }

        /* printf("val = %G finished newu = \n",*val); dprint(du,umin); */
        /* printf("umin = "); dprint(du,umin); */
    }
    else if (justrand == -1){
        // then start at each corner
        /* printf("corners\n"); */
        if (du == 1){
            nrand = 0;
            double distx = (ubu[0]-lbu[0])/4.0;
            ucurr[0] = lbu[0] + distx;
            res = c3opt_minimize(opt,ucurr,&valtemp);
            /* if (res < 0){ */
            /*     fprintf(stderr,"Optimization result %d\n",res); */
            /*     fprintf(stderr,"\t in bellman_optimal, du=1,justran=-1\n"); */

            /*     dprint(du,u); */
            /*     fprintf(stderr,"\n\n\n Rerunning with full output\n"); */

            /*     ucurr[0] = lbu[0] + distx; */
        
            /*     c3opt_set_verbose(opt,5); */
            /*     int res2 = c3opt_minimize(opt,ucurr,&valtemp); */
            /*     fprintf(stderr,"New optimization result %d\n",res2); */

            /*     diagnose_opt_failure(param,u,du); */
            /*     exit(1); */
            /* } */
            /* printf("val = %G, pt = ",valtemp);dprint(du,ucurr); */
            if (valtemp < minval){
                /* printf("\t updating current min!\n"); */
                minval = valtemp;
                memmove(umin,ucurr,du*sizeof(double));
            }
            ucurr[0] = ubu[0] - distx;
            c3opt_minimize(opt,ucurr,&valtemp);
            /* printf("val = %G, pt = ",valtemp);dprint(du,ucurr); */
            if (valtemp < minval){
                minval = valtemp;
                memmove(umin,ucurr,du*sizeof(double));
            }
        }
        else if (du == 2){
            nrand = 1; // 2
            npert = 1; // 1
            /* nrand = 2; */
            /* npert = 1; */
            double ut[2];
            for (size_t ii = 0; ii < 2; ii++){
                double distx = (ubu[0]-lbu[0])/4.0;
                for (size_t jj = 0; jj < 2; jj++){
                    double disty = (ubu[1]-lbu[1])/4.0;
                    if (ii == 0){
                        ucurr[0] = lbu[0]+distx;
                    }
                    else{
                        ucurr[0] = ubu[0]-distx;
                    }
                    if (jj == 0){
                        ucurr[1] = lbu[1]+disty;
                    }
                    else{
                        ucurr[1] = ubu[1]-disty;
                    }

                    c3opt_minimize(opt,ucurr,&valtemp);
                    /* printf("val = %G, pt = ",valtemp);dprint(du,ucurr); */
                    double dd = 0.0; // sufficient decrease
                    if (valtemp < (minval-dd)){
                        /* printf("\t updating current min!\n"); */
                        minval = valtemp;
                        memmove(umin,ucurr,du*sizeof(double));
                    }
                    for (size_t ll = 0; ll < npert; ll++){ // randomly perturb minimum
                        ut[0] = ucurr[0] + randu()*(distx*2.0) - distx;
                        ut[1] = ucurr[1] + randu()*(disty*2.0) - disty;
                        for (size_t zz = 0; zz < du; zz++){
                            if (ut[zz] > ubu[zz]){
                                ut[zz] = ubu[zz];
                            }
                            else if (ut[zz] < lbu[zz]){
                                ut[zz] = lbu[zz];
                            }
                        }
                        c3opt_minimize(opt,ut,&valtemp);
                        if (valtemp < (minval-dd)){
                            /* printf("\t updating current min!\n"); */
                            minval = valtemp;
                            memmove(umin,ut,du*sizeof(double));
                        }
                    }
                }
            }
        }
        else if (du == 3){
            nrand = 1; // 2
            npert = 1; // 1
            /* nrand = 2; */
            /* npert = 1; */
            double ut[3];
            for (size_t ii = 0; ii < 2; ii++){
                double distx = (ubu[0]-lbu[0])/4.0;
                for (size_t jj = 0; jj < 2; jj++){
                    double disty = (ubu[1]-lbu[1])/4.0;
                    for (size_t kk = 0; kk < 2; kk++){
                        double distz = (ubu[2]-lbu[2])/4.0;
                        if (ii == 0){
                            ucurr[0] = lbu[0]+distx;
                        }
                        else{
                            ucurr[0] = ubu[0]-distx;
                        }
                        if (jj == 0){
                            ucurr[1] = lbu[1]+disty;
                        }
                        else{
                            ucurr[1] = ubu[1]-disty;
                        }
                        if (kk == 0){
                            ucurr[2] = lbu[2]+distz;
                        }
                        else{
                            ucurr[2] = ubu[2]-distz;
                        }
                        c3opt_minimize(opt,ucurr,&valtemp);
                        /* printf("val = %G, pt = ",valtemp);dprint(du,ucurr); */
                        double dd = 0.0; // sufficient decrease
                        if (valtemp < (minval-dd)){
                            /* printf("\t updating current min!\n"); */
                            minval = valtemp;
                            memmove(umin,ucurr,du*sizeof(double));
                        }
                        for (size_t ll = 0; ll < npert; ll++){ // randomly perturb minimum
                            ut[0] = ucurr[0] + randu()*(distx*2.0) - distx;
                            ut[1] = ucurr[1] + randu()*(disty*2.0) - disty;
                            ut[2] = ucurr[2] + randu()*(distz*2.0) - distz;
                            for (size_t zz = 0; zz < du; zz++){
                                if (ut[zz] > ubu[zz]){
                                    ut[zz] = ubu[zz];
                                }
                                else if (ut[zz] < lbu[zz]){
                                    ut[zz] = lbu[zz];
                                }
                            }
                            c3opt_minimize(opt,ut,&valtemp);
                            if (valtemp < (minval-dd)){
                                /* printf("\t updating current min!\n"); */
                                minval = valtemp;
                                memmove(umin,ut,du*sizeof(double));
                            }
                        }
                    }
                }
            }
        }
        else if (du == 4){
            nrand = 1; // 2
            npert = 1; // 1
            /* nrand = 2; */
            /* npert = 1; */
            double ut[4];
            for (size_t ii = 0; ii < 2; ii++){
                double distx = (ubu[0]-lbu[0])/4.0;
                for (size_t jj = 0; jj < 2; jj++){
                    double disty = (ubu[1]-lbu[1])/4.0;
                    for (size_t kk = 0; kk < 2; kk++){
                        double distz = (ubu[2]-lbu[2])/4.0;
                        for (size_t aa = 0; aa < 2; aa++){
                            double dista = (ubu[3]-lbu[3])/4.0;
                            if (ii == 0){
                                ucurr[0] = lbu[0]+distx;
                            }
                            else{
                                ucurr[0] = ubu[0]-distx;
                            }
                            if (jj == 0){
                                ucurr[1] = lbu[1]+disty;
                            }
                            else{
                                ucurr[1] = ubu[1]-disty;
                            }
                            if (kk == 0){
                                ucurr[2] = lbu[2]+distz;
                            }
                            else{
                                ucurr[2] = ubu[2]-distz;
                            }
                            if (aa == 0){
                                ucurr[3] = lbu[3]+dista;
                            }
                            else{
                                ucurr[3] = ubu[3]-dista;
                            }
                            c3opt_minimize(opt,ucurr,&valtemp);
                            /* printf("val = %G, pt = ",valtemp);dprint(du,ucurr); */
                            double dd = 0.0; // sufficient decrease
                            if (valtemp < (minval-dd)){
                                /* printf("\t updating current min!\n"); */
                                minval = valtemp;
                                memmove(umin,ucurr,du*sizeof(double));
                            }
                            for (size_t ll = 0; ll < npert; ll++){ // randomly perturb minimum
                                ut[0] = ucurr[0] + randu()*(distx*2.0) - distx;
                                ut[1] = ucurr[1] + randu()*(disty*2.0) - disty;
                                ut[2] = ucurr[2] + randu()*(distz*2.0) - distz;
                                ut[3] = ucurr[3] + randu()*(dista*2.0) - dista;
                                for (size_t zz = 0; zz < du; zz++){
                                    if (ut[zz] > ubu[zz]){
                                        ut[zz] = ubu[zz];
                                    }
                                    else if (ut[zz] < lbu[zz]){
                                        ut[zz] = lbu[zz];
                                    }
                                }
                                c3opt_minimize(opt,ut,&valtemp);
                                if (valtemp < (minval-dd)){
                                    /* printf("\t updating current min!\n"); */
                                    minval = valtemp;
                                    memmove(umin,ut,du*sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (du == 5){
            nrand = 1; // 2
            npert = 1; // 1
            /* nrand = 2; */
            /* npert = 1; */
            double ut[5];
            for (size_t ii = 0; ii < 2; ii++){
                double distx = (ubu[0]-lbu[0])/4.0;
                for (size_t jj = 0; jj < 2; jj++){
                    double disty = (ubu[1]-lbu[1])/4.0;
                    for (size_t kk = 0; kk < 2; kk++){
                        double distz = (ubu[2]-lbu[2])/4.0;
                        for (size_t aa = 0; aa < 2; aa++){
                            double dista = (ubu[3]-lbu[3])/4.0;
                            for (size_t bb = 0; bb < 2; bb++){
                                double distb = (ubu[4]-lbu[4])/4.0;
                                if (ii == 0){
                                    ucurr[0] = lbu[0]+distx;
                                }
                                else{
                                    ucurr[0] = ubu[0]-distx;
                                }
                                if (jj == 0){
                                    ucurr[1] = lbu[1]+disty;
                                }
                                else{
                                    ucurr[1] = ubu[1]-disty;
                                }
                                if (kk == 0){
                                    ucurr[2] = lbu[2]+distz;
                                }
                                else{
                                    ucurr[2] = ubu[2]-distz;
                                }
                                if (aa == 0){
                                    ucurr[3] = lbu[3]+dista;
                                }
                                else{
                                    ucurr[3] = ubu[3]-dista;
                                }
                                if (bb == 0){
                                    ucurr[4] = lbu[4]+distb;
                                }
                                else{
                                    ucurr[4] = ubu[4]-distb;
                                }
                                c3opt_minimize(opt,ucurr,&valtemp);
                                /* printf("val = %G, pt = ",valtemp);dprint(du,ucurr); */
                                double dd = 0.0; // sufficient decrease
                                if (valtemp < (minval-dd)){
                                    /* printf("\t updating current min!\n"); */
                                    minval = valtemp;
                                    memmove(umin,ucurr,du*sizeof(double));
                                }
                                for (size_t ll = 0; ll < npert; ll++){ // randomly perturb minimum
                                    ut[0] = ucurr[0] + randu()*(distx*2.0) - distx;
                                    ut[1] = ucurr[1] + randu()*(disty*2.0) - disty;
                                    ut[2] = ucurr[2] + randu()*(distz*2.0) - distz;
                                    ut[3] = ucurr[3] + randu()*(dista*2.0) - dista;
                                    ut[4] = ucurr[4] + randu()*(distb*2.0) - distb;
                                    for (size_t zz = 0; zz < du; zz++){
                                        if (ut[zz] > ubu[zz]){
                                            ut[zz] = ubu[zz];
                                        }
                                        else if (ut[zz] < lbu[zz]){
                                            ut[zz] = lbu[zz];
                                        }
                                    }
                                    c3opt_minimize(opt,ut,&valtemp);
                                    if (valtemp < (minval-dd)){
                                        /* printf("\t updating current min!\n"); */
                                        minval = valtemp;
                                        memmove(umin,ut,du*sizeof(double));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (du == 6){
            nrand = 1; // 2
            npert = 1; // 1
            /* nrand = 2; */
            /* npert = 1; */
            double ut[6];
            for (size_t ii = 0; ii < 2; ii++){
                double distx = (ubu[0]-lbu[0])/4.0;
                for (size_t jj = 0; jj < 2; jj++){
                    double disty = (ubu[1]-lbu[1])/4.0;
                    for (size_t kk = 0; kk < 2; kk++){
                        double distz = (ubu[2]-lbu[2])/4.0;
                        for (size_t aa = 0; aa < 2; aa++){
                            double dista = (ubu[3]-lbu[3])/4.0;
                            for (size_t bb = 0; bb < 2; bb++){
                                double distb = (ubu[4]-lbu[4])/4.0;
                                for (size_t cc = 0; cc < 2; cc++){
                                    double distc = (ubu[5]-lbu[5])/4.0;
                                    if (ii == 0){
                                        ucurr[0] = lbu[0]+distx;
                                    }
                                    else{
                                        ucurr[0] = ubu[0]-distx;
                                    }
                                    if (jj == 0){
                                        ucurr[1] = lbu[1]+disty;
                                    }
                                    else{
                                        ucurr[1] = ubu[1]-disty;
                                    }
                                    if (kk == 0){
                                        ucurr[2] = lbu[2]+distz;
                                    }
                                    else{
                                        ucurr[2] = ubu[2]-distz;
                                    }
                                    if (aa == 0){
                                        ucurr[3] = lbu[3]+dista;
                                    }
                                    else{
                                        ucurr[3] = ubu[3]-dista;
                                    }
                                    if (bb == 0){
                                        ucurr[4] = lbu[4]+distb;
                                    }
                                    else{
                                        ucurr[4] = ubu[4]-distb;
                                    }
                                    if (cc == 0){
                                        ucurr[5] = lbu[5]+distc;
                                    }
                                    else{
                                        ucurr[5] = ubu[5]-distc;
                                    }
                                    c3opt_minimize(opt,ucurr,&valtemp);
                                    /* printf("val = %G, pt = ",valtemp);dprint(du,ucurr); */
                                    double dd = 0.0; // sufficient decrease
                                    if (valtemp < (minval-dd)){
                                        /* printf("\t updating current min!\n"); */
                                        minval = valtemp;
                                        memmove(umin,ucurr,du*sizeof(double));
                                    }
                                    for (size_t ll = 0; ll < npert; ll++){ // randomly perturb minimum
                                        ut[0] = ucurr[0] + randu()*(distx*2.0) - distx;
                                        ut[1] = ucurr[1] + randu()*(disty*2.0) - disty;
                                        ut[2] = ucurr[2] + randu()*(distz*2.0) - distz;
                                        ut[3] = ucurr[3] + randu()*(dista*2.0) - dista;
                                        ut[4] = ucurr[4] + randu()*(distb*2.0) - distb;
                                        ut[5] = ucurr[5] + randu()*(distc*2.0) - distc;
                                        for (size_t zz = 0; zz < du; zz++){
                                            if (ut[zz] > ubu[zz]){
                                                ut[zz] = ubu[zz];
                                            }
                                            else if (ut[zz] < lbu[zz]){
                                                ut[zz] = lbu[zz];
                                            }
                                        }
                                        c3opt_minimize(opt,ut,&valtemp);
                                        if (valtemp < (minval-dd)){
                                            /* printf("\t updating current min!\n"); */
                                            minval = valtemp;
                                            memmove(umin,ut,du*sizeof(double));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (du == 7){
            nrand = 1; // 2
            npert = 1; // 1
            /* nrand = 2; */
            /* npert = 1; */
            double ut[7];
            for (size_t ii = 0; ii < 2; ii++){
                double distx = (ubu[0]-lbu[0])/4.0;
                for (size_t jj = 0; jj < 2; jj++){
                    double disty = (ubu[1]-lbu[1])/4.0;
                    for (size_t kk = 0; kk < 2; kk++){
                        double distz = (ubu[2]-lbu[2])/4.0;
                        for (size_t aa = 0; aa < 2; aa++){
                            double dista = (ubu[3]-lbu[3])/4.0;
                            for (size_t bb = 0; bb < 2; bb++){
                                double distb = (ubu[4]-lbu[4])/4.0;
                                for (size_t cc = 0; cc < 2; cc++){
                                    double distc = (ubu[5]-lbu[5])/4.0;
                                    for (size_t dd = 0; dd < 2; dd++){
                                        double distd = (ubu[6]-lbu[6])/4.0;
                                        if (ii == 0){
                                            ucurr[0] = lbu[0]+distx;
                                        }
                                        else{
                                            ucurr[0] = ubu[0]-distx;
                                        }
                                        if (jj == 0){
                                            ucurr[1] = lbu[1]+disty;
                                        }
                                        else{
                                            ucurr[1] = ubu[1]-disty;
                                        }
                                        if (kk == 0){
                                            ucurr[2] = lbu[2]+distz;
                                        }
                                        else{
                                            ucurr[2] = ubu[2]-distz;
                                        }
                                        if (aa == 0){
                                            ucurr[3] = lbu[3]+dista;
                                        }
                                        else{
                                            ucurr[3] = ubu[3]-dista;
                                        }
                                        if (bb == 0){
                                            ucurr[4] = lbu[4]+distb;
                                        }
                                        else{
                                            ucurr[4] = ubu[4]-distb;
                                        }
                                        if (cc == 0){
                                            ucurr[5] = lbu[5]+distc;
                                        }
                                        else{
                                            ucurr[5] = ubu[5]-distc;
                                        }
                                        if (dd == 0){
                                            ucurr[6] = lbu[6]+distc;
                                        }
                                        else{
                                            ucurr[6] = ubu[6]-distc;
                                        }
                                        c3opt_minimize(opt,ucurr,&valtemp);
                                        /* printf("val = %G, pt = ",valtemp);dprint(du,ucurr); */
                                        double dec = 0.0; // sufficient decrease
                                        if (valtemp < (minval-dec)){
                                            /* printf("\t updating current min!\n"); */
                                            minval = valtemp;
                                            memmove(umin,ucurr,du*sizeof(double));
                                        }
                                        for (size_t ll = 0; ll < npert; ll++){ // randomly perturb minimum
                                            ut[0] = ucurr[0] + randu()*(distx*2.0) - distx;
                                            ut[1] = ucurr[1] + randu()*(disty*2.0) - disty;
                                            ut[2] = ucurr[2] + randu()*(distz*2.0) - distz;
                                            ut[3] = ucurr[3] + randu()*(dista*2.0) - dista;
                                            ut[4] = ucurr[4] + randu()*(distb*2.0) - distb;
                                            ut[5] = ucurr[5] + randu()*(distc*2.0) - distc;
                                            ut[6] = ucurr[6] + randu()*(distd*2.0) - distd;
                                            for (size_t zz = 0; zz < du; zz++){
                                                if (ut[zz] > ubu[zz]){
                                                    ut[zz] = ubu[zz];
                                                }
                                                else if (ut[zz] < lbu[zz]){
                                                    ut[zz] = lbu[zz];
                                                }
                                            }
                                            c3opt_minimize(opt,ut,&valtemp);
                                            if (valtemp < (minval-dec)){
                                                /* printf("\t updating current min!\n"); */
                                                minval = valtemp;
                                                memmove(umin,ut,du*sizeof(double));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    
        for (size_t jj = 0; jj < nrand; jj++){
            for (size_t kk = 0; kk < du; kk++){
                ucurr[kk] = randu()*(ubu[kk]-lbu[kk]) + lbu[kk];
            }
            c3opt_minimize(opt,ucurr,&valtemp);
            if (valtemp < minval){
                /* printf("\t updating current min!\n"); */
                minval = valtemp;
                memmove(umin,ucurr,du*sizeof(double));
            }
        }
    }
    
    memmove(u,umin,du*sizeof(double));
    *val = minval;

    c3opt_free(opt); opt = NULL;
    free(umin); umin = NULL;
    free(ucurr); ucurr = NULL;
    return 0;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

struct VIparam
{
    struct ControlParams * cp;
    struct ValueF * vf;

    size_t nstate_evals;
    size_t nnode_evals;
    double convergence;
    double time_in_loop;

    size_t oniter;
};

struct VIparam * vi_param_create(double convergence)
{
    struct VIparam * vi = malloc(sizeof(struct VIparam));
    assert (vi != NULL);
    vi->convergence = convergence;
    
    vi->nstate_evals = 0; // counts in terms of fibers
    vi->nnode_evals  = 0; // I think this is redundant
    vi->cp = NULL;
    vi->vf = NULL;
    vi->time_in_loop = 0;
    
    return vi;
}

void vi_param_destroy(struct VIparam * vi)
{
    if (vi != NULL){
        free(vi); vi = NULL;
    }
}

void vi_param_add_cp(struct VIparam * vi, struct ControlParams * cp)
{
    assert (vi != NULL);
    vi->cp = cp;
}

void vi_param_add_value(struct VIparam * vi, struct ValueF * vf)
{
    assert (vi != NULL);
    vi->vf = vf;

    vi->nstate_evals = 0;
    vi->nnode_evals = 0;
}

/* int bellman_vi(size_t N, const double * x, double * out, void * arg) */
/* { */
/*     struct VIparam * param = arg; */
/*     struct ControlParams * cp = param->cp; */
/*     struct ValueF * vf = param->vf; */
/*     assert (cp != NULL); */
/*     assert (vf != NULL); */
/*     struct MCAparam * mca = cp->mca; */
/*     struct DPparam * dp = cp->dp; */

/*     size_t dx = mca->dx; */
/*     size_t du = mca->du; */
/*     size_t * ngrid = mca->ngrid; */
/*     double ** xgrid = mca->xgrid; */
    

/*     struct HTable * htable = workspace_get_vi_htable(cp->work); */
/*     struct Boundary * bound = dp->bound; */

/*     int * absorbed = workspace_get_absorbed(cp->work,0); */
/*     double * costs = workspace_get_costs(cp->work,0); */

/*     /\* printf("get neighbor cost\n"); *\/ */
/*     size_t * fi = calloc_size_t(dx); */

/*     size_t dim_vary; */
/*     int res = mca_get_neighbor_costs(dx,N,x, bound,vf,ngrid,xgrid, */
/*                                      fi,&dim_vary, absorbed,costs); */

/*     assert (res == 0); */
/*     size_t * ind_to_serialize = workspace_get_ind_to_serialize(cp->work); */
/*     for (size_t ii = 0; ii < dx; ii++){ */
/*         if (ii != dim_vary){ */
/*             ind_to_serialize[ii] = fi[ii]; */
/*         } */
/*     } */
/*     ind_to_serialize[dx] = dim_vary; */
/*     ind_to_serialize[dx+1] = workspace_get_vi_iter(cp->work); */
    
/*     char * key = size_t_a_to_char(ind_to_serialize,dx+2); */
/*     size_t nbytes = 0; */
/*     double * out_stored = htable_get_element(htable,key,&nbytes); */
/*     if (out_stored == NULL){ */
/*         double time = 0.0; */
/*         control_params_add_time_and_states(cp,time,N,x); */
/*         param->nstate_evals = param->nstate_evals + N; */

/*         size_t * absorbed_no = workspace_get_absorbed_no(cp->work); */
/*         size_t * absorbed_yes = workspace_get_absorbed_yes(cp->work); */
/*         size_t ano_ind = 0; */
/*         size_t ayes_ind = 0; */
/*         for (size_t ii = 0; ii < N; ii++){ */
/*             if (absorbed[ii] == 0){ */
/*                 absorbed_no[ano_ind] = ii; */
/*                 ano_ind++; */
/*             } */
/*             else{ */
/*                 absorbed_yes[ayes_ind] = ii; */
/*                 ayes_ind++;; */
/*             } */
/*         } */

/*         // this loop is very fast because only deals with absorbed things */
/*         for (size_t ii = 0; ii < ayes_ind; ii++){ */
/*             struct Memory mem; */
/*             mem.shared = cp; */
/*             mem.private = absorbed_yes[ii]; */
/*             double * u = workspace_get_u(cp->work,absorbed_yes[ii]); */
/*             bellman_optimal(du,u,out+absorbed_yes[ii],&mem); */
/*         } */
        
/*         /\* printf("\n\n\n\n"); *\/ */
/*         double t_start = omp_get_wtime(); */
/*         #pragma omp parallel for schedule(guided) */
/*         for (size_t ii = 0; ii < ano_ind; ii++){ */
/*             /\* printf("Hello from thread %d/%d\n",omp_get_thread_num(),omp_get_num_threads()); *\/ */
/*             struct Memory mem; */
/*             mem.shared = cp; */
/*             mem.private = absorbed_no[ii]; */
            
/*             double * u = workspace_get_u(cp->work,absorbed_no[ii]); */
/*             double val; */
/*             /\* double t_start_2 = omp_get_wtime(); *\/ */
/*             int res2 = bellman_optimal(du,u,&val,&mem); */
/*             /\* double t_end_2 = omp_get_wtime(); *\/ */
/*             /\* printf("Timing: %zu ,%d, %G\n",absorbed_no[ii],absorbed[absorbed_no[ii]],t_end_2-t_start_2); *\/ */
/*             /\* printf("\t x = (%G,%G)\n\n\n",x[ absorbed_no[ii]*dx], x[absorbed_no[ii]*dx + 1] ); *\/ */
/*             assert (res2 == 0); */
/*             out[absorbed_no[ii]] = val; */
/*         } */
/*         double t_end = omp_get_wtime(); */
/*         double run_time = (double)(t_end - t_start); /\* / CLOCKS_PER_SEC; *\/ */
/*         param->time_in_loop += run_time; */
/*         /\* printf("%G\n",run_time); *\/ */

/*         param->nnode_evals += N; */
/*         htable_add_element(htable,key,out,N * sizeof(double)); */
/*     } */
/*     else{ */
/*         for (size_t ii = 0; ii < N; ii++){ */
/*             out[ii] = out_stored[ii]; */
/*         } */
/*         free(key); key = NULL; */
/*     } */

/*     free(fi); fi = NULL; */


/*     return 0; */
    
/* } */


int bellman_vi(size_t N, const double * x, double * out, void * arg)
{
    struct VIparam * param = arg;
    struct ControlParams * cp = param->cp;
    struct ValueF * vf = param->vf;
    assert (cp != NULL);
    assert (vf != NULL);
    struct MCAparam * mca = cp->mca;
    struct DPparam * dp = cp->dp;

    size_t dx = mca->dx;
    size_t du = mca->du;
    size_t * ngrid = mca->ngrid;
    double ** xgrid = mca->xgrid;
    

    struct HTable * htable = workspace_get_vi_htable(cp->work);
    struct Boundary * bound = dp->bound;

    int * absorbed = workspace_get_absorbed(cp->work,0);
    double * costs = workspace_get_costs(cp->work,0);

    /* printf("get neighbor cost\n"); */
    size_t * fi = calloc_size_t(dx);

    size_t dim_vary;
    int res = mca_get_neighbor_costs(dx,N,x, bound,vf,ngrid,xgrid,
                                     fi,&dim_vary, absorbed,costs);

    double time = 0.0;
    control_params_add_time_and_states(cp,time,N,x);
    
    size_t nrun = 0;
    size_t * ind_run = calloc_size_t(N);
    char ** saved_keys = workspace_get_saved_keys(cp->work);/* malloc(N * sizeof(char*)); */
    assert (saved_keys != NULL);

    assert (res == 0);
    size_t * ind_to_serialize = workspace_get_ind_to_serialize(cp->work);
    for (size_t ii = 0; ii < dx; ii++){
        ind_to_serialize[ii] = fi[ii];
    }
    ind_to_serialize[dx] = 0;/* dim_vary; */
    ind_to_serialize[dx+1] = workspace_get_vi_iter(cp->work);

    /* printf("\n"); */
    for (size_t ii = 0; ii < N; ii++){
        ind_to_serialize[dim_vary] = ii;
        
        size_t_a_to_char(ind_to_serialize,dx+2,saved_keys[ii]);
        /* printf("serializing: "); iprint_sz(dx+2,ind_to_serialize); */
        /* printf("x = "); dprint(dx,x+ii*dx); */
        /* printf("key = %s\n",saved_keys[ii]); */
        size_t nbytes = 0;
        double * out_stored = htable_get_element(htable,saved_keys[ii],&nbytes);

        /* printf("nbytes = %zu,%zu\n",nbytes,sizeof(double)); */
        if (out_stored != NULL){
            out[ii] = out_stored[0];

            /* struct Memory mem;  */
            /* mem.shared = cp; */
            /* mem.private = ii; */
            /* double val; */
            /* double * u = workspace_get_u(cp->work,ii); */
            /* bellman_optimal(du,u,&val,&mem); */
            /* double diff = val - out[ii]; */
            /* if (fabs(diff) > 1e-14){ */
            /*     printf("out = %3.15G\n",out[ii]); */
            /*     printf("val = %3.15G\n",val); */
            /*     printf("diff = %G\n",val-out[ii]); */
            /*     printf("ii = %zu\n",ii); */
            /*     dprint(dx,x + ii * dx); */
            /*     exit(1); */
            /* } */
        }
        else if (absorbed[ii] == 0){ // not absorbed
            ind_run[nrun] = ii;
            nrun++;
            param->nstate_evals++; // these are the ones which I will run below
            param->nnode_evals++;
        }
        else{
            struct Memory mem; 
            mem.shared = cp;
            mem.private = ii;
            double * u = workspace_get_u(cp->work,ii);
            bellman_optimal(du,u,out+ii,&mem);
            htable_add_element(htable,saved_keys[ii],out+ii,1);
            /* printf("out[%zu] = %G\n",ii,out[ii]); */
            param->nstate_evals++; // ran it
            param->nnode_evals++;
        }
    }

    #ifdef _OPENMP
    double t_start = omp_get_wtime();
    #pragma omp parallel for schedule(guided)
    #endif
    for (size_t ii = 0; ii < nrun; ii++){
        struct Memory mem;
        mem.shared = cp;
        mem.private = ind_run[ii];
            
        double * u = workspace_get_u(cp->work,ii);
        int res2 = bellman_optimal(du,u,out + ind_run[ii],&mem);
        /* htable_add_element(htable,saved_keys[ind_run[ii]],&out[ind_run[ii]],1); */
        assert (res2 == 0);
        /* out[ind_run[ii]] = val; */
    }
    #ifdef _OPENMP
    double t_end = omp_get_wtime();
    double run_time = (double)(t_end - t_start); /* / CLOCKS_PER_SEC; */
    param->time_in_loop += run_time;
    #endif


    // doing this not in parallel because not thread safe
    for (size_t ii = 0; ii < nrun; ii++){
        char * key = saved_keys[ind_run[ii]];
        /* printf("key = %s, out = %3.5G\n",key,out[ind_run[ii]]); */
        htable_add_element(htable,key,&out[ind_run[ii]],1);
    }

    free(fi); fi = NULL;
    free(ind_run); ind_run = NULL;
    /* free(saved_keys); saved_keys = NULL; */
    return 0;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
struct PIparam
{
    struct ControlParams * cp;
    struct ValueF * vf_iteration;
    /* struct HTable * htable_iter; // store fiber evaluations for fixed iteration */

    struct ValueF * vf_policy;
    /* struct HTable * htable; // store probability transitions, stage costs, etc associated with policy */

    size_t npol_evals; // evaluation of policies
    size_t niter_evals; // evaluation of fibers for one step of PI
    size_t niter_node_evals; // evaluation of nodes for one step of VI
    double convergence;
};

struct PIparam * pi_param_create(double convergence, struct ValueF * policy)
{
    struct PIparam * poli = malloc(sizeof(struct PIparam));
    assert (poli != NULL);
    poli->convergence = convergence;
    poli->cp = NULL;
    
    poli->vf_policy = policy;
    /* poli->htable = htable_create(100000); */
    poli->npol_evals = 0;
        
    poli->vf_iteration = NULL;
    
    /* poli->htable_iter = NULL; */
    poli->niter_evals = 0;
    poli->niter_node_evals = 0;

    return poli;
}

void pi_param_destroy(struct PIparam * poli)
{
    if (poli != NULL){

        /* htable_destroy(poli->htable); poli->htable = NULL; */
        /* htable_destroy(poli->htable_iter); poli->htable_iter = NULL; */
        free(poli); poli = NULL;
    }
}

void pi_param_add_cp(struct PIparam * poli, struct ControlParams * cp)
{
    assert (poli != NULL);
    poli->cp = cp;
}

void pi_param_add_value(struct PIparam * poli, struct ValueF * vf)
{
    assert (poli != NULL);
    poli->vf_iteration = vf;

    /* htable_destroy(poli->htable_iter);  */
    /* poli->htable_iter = NULL; */
    /* poli->htable_iter = htable_create(1000000); */
    poli->niter_evals = 0;
    poli->niter_node_evals = 0;
}

/* int bellman_pi(size_t N, const double * x, double * out, void * arg) */
/* { */

/*     /\* printf("started pi\n"); *\/ */
/*     struct PIparam * param = arg; */
/*     assert (param->cp != NULL); */
/*     struct ControlParams * cp = param->cp; */
/*     struct MCAparam * mca = cp->mca; */
/*     assert (mca != NULL); */
/*     struct DPparam * dp = cp->dp; */
/*     assert (dp != NULL); */
/*     size_t dx = mca->dx; */
/*     size_t du = mca->du; */
/*     size_t * ngrid = mca->ngrid; */
/*     double ** xgrid = mca->xgrid; */
    
/*     struct ValueF * vf_policy    = param->vf_policy; */
/*     struct ValueF * vf_iteration = param->vf_iteration; */
/*     assert (vf_policy != NULL); */
/*     assert (vf_iteration != NULL); */
/*     struct Boundary * bound = dp->bound; */
/*     assert (bound != NULL); */

/*     /\* struct HTable * htable = param->htable; *\/ */
/*     /\* struct HTable * htable_iter = param->htable_iter; *\/ */

/*     struct HTable * htable = workspace_get_pi_prob_htable(cp->work); */
/*     struct HTable * htable_iter = workspace_get_pi_htable(cp->work); */

/*     assert (htable != NULL); */
/*     assert (htable_iter != NULL); */
/*     // handle the policy */

/*     double time = 0.0; */
/*     int * absorbed_pol = workspace_get_absorbed(cp->work,0); */
/*     double * costs_pol = workspace_get_costs(cp->work,0); */

/*     size_t * fi = calloc_size_t(dx); */
/*     size_t dim_vary; */
/*     int res = mca_get_neighbor_costs(dx,N,x,bound,vf_policy,ngrid,xgrid, */
/*                                      fi,&dim_vary,absorbed_pol,costs_pol); */

/*     assert (res == 0); */
/*     assert (dim_vary < dx); */

/*     size_t * ind_to_serialize = workspace_get_ind_to_serialize(cp->work); */
/*     for (size_t ii = 0; ii < dx; ii++){ */
/*         assert (fi[ii] < ngrid[ii]); */
/*         if (ii != dim_vary){ */
/*             ind_to_serialize[ii] = fi[ii]; */
/*         } */
/*     } */
/*     ind_to_serialize[dx] = dim_vary; */
/*     ind_to_serialize[dx+1] = workspace_get_pi_iter(cp->work); // stores values for probabilities for current policy */
/*     ind_to_serialize[dx+2] = workspace_get_pi_subiter(cp->work); // stores value for the current cost function (changes faster than policy) */
/*     free(fi); fi = NULL; */
/*     char * key1 = size_t_a_to_char(ind_to_serialize,dx+3); */
/*     char * key2 = size_t_a_to_char(ind_to_serialize,dx+2); // NOTE DIFFERENCE, only uses major iter not subiter */
/*     char * key3 = size_t_a_to_char(ind_to_serialize,dx+3); */

/*     // see if already computed the output for this function */
/*     size_t nbytes1 = 0; */
/*     double * out_stored = htable_get_element(htable_iter,key1,&nbytes1); */
/*     if (out_stored == NULL){ */
/*         /\* printf("not stored\n"); *\/ */
/*         size_t nbytes = 0; */
/*         // see if already computed transition probabilities for this fiber */
/*         double * probs = htable_get_element(htable,key2,&nbytes); */
/*         /\* probs = NULL; *\/ */
/*         /\* printf("got probs\n"); *\/ */
/*         int probs_alloc = 0; */
/*         /\* probs = NULL; *\/ */
/*         if (probs == NULL){ */
/*             /\* printf("probs is null\n");; *\/ */
/*             size_t nprobs = N*(2*dx+1) + 2*N; // last N are stage costs and N before that are dts */
/*             probs = calloc_double(nprobs); // last N are dts */
/*             probs_alloc = 1; */

/*             control_params_add_time_and_states(cp,time,N,x); */

/*             #pragma omp parallel for schedule(guided) */
/*             /\* #pragma omp parallel for schedule(static) *\/ */
/*             for (size_t ii = 0; ii < N; ii++){ */

/*                 struct Memory mem; */
/*                 mem.shared = cp; */
/*                 mem.private = ii; */
                

/*                 double * u = workspace_get_u(cp->work,ii); */

/*                 double val; */
/*                 int res2 = bellman_optimal(du,u,&val,&mem); */

                
/*                 double * drift = workspace_get_drift(cp->work, ii); */
/*                 double * diff  = workspace_get_diff (cp->work, ii); */


/*                 // build transition probabilities at optimal */
/*                 assert (res2 == 0); */
/*                 res2 = drift_eval(dp->dyn_drift,time,x+ii*dx,u,drift,NULL); */
/*                 assert (res2 == 0); */
/*                 res2 = diff_eval(dp->dyn_diff,time,x+ii*dx,u,diff,NULL); */
/*                 /\* res = transition_assemble(dx,du,cp->dw,mca->hmin,mca->hvec, *\/ */
/*                 /\*                           drift,NULL,diff,NULL, *\/ */
/*                 /\*                           probs + ii*(2*dx+1), *\/ */
/*                 /\*                           NULL,probs + N*(2*dx+1)+ii, // dt *\/ */
/*                 /\*                           NULL,NULL); *\/ */
/*                 res = transition_assemble(dx,du,cp->dw,mca->h2,mca->t, */
/*                                           drift,NULL,diff,NULL, */
/*                                           probs + ii*(2*dx+1), */
/*                                           NULL,probs + N*(2*dx+1)+ii, // dt */
/*                                           NULL,NULL); */

/*                 res2 = dp->stagecost(time,x+ii*dx,u,probs + N*(2*dx+1)+N+ii,NULL); */
/*                 assert (res2 == 0); */


/*             } */
/*             param->npol_evals += N; */
/*             htable_add_element(htable,key2,probs,nprobs * sizeof(double)); */

/*         } */
/*         else{ */
/*             /\* printf("probs is null"); *\/ */
/*             free(key2); key2 = NULL; */
/*         } */

/*         // need to compute the control somewhere */
/*         int * absorbed = calloc_int(N); */
/*         double * costs = calloc_double(N*(2*dx+1)); */
/*         fi = calloc_size_t(dx); */
/*         /\* printf("get other mca\n"); *\/ */
/*         res = mca_get_neighbor_costs(dx,N,x, */
/*                                      bound,vf_iteration,ngrid,xgrid, */
/*                                      fi,&dim_vary, */
/*                                      absorbed,costs); */
/*         /\* printf("got other \n"); iprint(N,absorbed); *\/ */
        
/*         assert (res == 0); */


/*         param->niter_evals = param->niter_evals + N; */
/*         param->niter_node_evals = param->niter_node_evals + N; */
/*         /\* printf("N = %zu\n",N); *\/ */
/*         for (size_t ii = 0; ii < N; ii++){ */

            
/*             // get the optimal control */
/*             // iterate with the optimal control */
/*             if (absorbed[ii] == 0){ */
/*                 out[ii] = bellmanrhs(dx,du, */
/*                                      *(probs + N*(2*dx+1)+N+ii), // stagecost */
/*                                      NULL, */
/*                                      dp->discount, */
/*                                      probs+ii*(2*dx+1),NULL, */
/*                                      *(probs+N*(2*dx+1)+ii),NULL, // dt */
/*                                      costs+ii*(2*dx+1),NULL); */
/*                 /\* printf("out[%zu] = %G\n",ii,out[ii]); *\/ */
             
/*             } */
/*             else if (absorbed[ii] == 1){ // absorbed cost */
/*                 /\* printf("boundary[%zu]\n",ii); *\/ */
/*                 /\* printf("\t"); dprint(dx,x+ii*dx); *\/ */
/*                 res = dp->boundcost(time,x+ii*dx,out+ii); */
/*                 assert (res == 0); */
/*             } */
/*             else if (absorbed[ii] == -1){ */
/*                 /\* printf("weird [%zu]\n",ii); *\/ */
/*                 /\* printf("\t"); dprint(dx,x+ii*dx); *\/ */
/*                 /\* assert (1 == 0); *\/ */
/*                 res = dp->obscost(x+ii*dx,out+ii); */
/*             } */
/*             else{ */
/*                 printf("ii = %zu, N = %zu, dx = %zu\n",ii,N,dx); */
/*                 printf("x = "); dprint(dx,x+ii*dx); */
/*                 printf("absorbed = "); iprint(N,absorbed); */
/*                 fprintf(stderr, "Unrecognized aborbed condition %d\n",absorbed[ii]); */
/*                 exit(1); */
/*             } */

/*         } */

/*         /\* printf("free stuff\n"); *\/ */
/*         if (probs_alloc == 1){ */
/*             free(probs); probs = NULL; */
/*         } */
/*         free(fi); fi = NULL; */

/*         htable_add_element(htable_iter,key3,out,N * sizeof(double)); */
/*     } */
/*     else{ */
/*         /\* printf("got saved\n"); *\/ */
/*         for (size_t ii = 0; ii < N; ii++){ */
/*             out[ii] = out_stored[ii]; */
/*         } */
/*         /\* dprint(N,out); *\/ */
/*         free(key3); key3 = NULL; */
/*         free(key2); key2 = NULL; */
/*     } */

/*     free(key1); key1 = NULL; */
/*     /\* free(key2); key2 = NULL; *\/ */

/*     return 0; */
    
/* } */

int bellman_pi(size_t N, const double * x, double * out, void * arg)
{

    /* printf("started pi\n"); */
    struct PIparam * param = arg;
    assert (param->cp != NULL);
    struct ControlParams * cp = param->cp;
    struct MCAparam * mca = cp->mca;
    assert (mca != NULL);
    struct DPparam * dp = cp->dp;
    assert (dp != NULL);
    size_t dx = mca->dx;
    size_t du = mca->du;
    size_t * ngrid = mca->ngrid;
    double ** xgrid = mca->xgrid;
    
    struct ValueF * vf_policy    = param->vf_policy;
    struct ValueF * vf_iteration = param->vf_iteration;
    assert (vf_policy != NULL);
    assert (vf_iteration != NULL);
    struct Boundary * bound = dp->bound;
    assert (bound != NULL);

    /* struct HTable * htable = param->htable; */
    /* struct HTable * htable_iter = param->htable_iter; */

    struct HTable * htable = workspace_get_pi_prob_htable(cp->work);
    struct HTable * htable_iter = workspace_get_pi_htable(cp->work);

    assert (htable != NULL);
    assert (htable_iter != NULL);
    // handle the policy

    double time = 0.0;
    control_params_add_time_and_states(cp,time,N,x);
    int * absorbed_pol = workspace_get_absorbed(cp->work,0);
    double * costs_pol = workspace_get_costs(cp->work,0);

    size_t * fi = calloc_size_t(dx);
    size_t dim_vary;
    int res = mca_get_neighbor_costs(dx,N,x,bound,vf_policy,ngrid,xgrid,
                                     fi,&dim_vary,absorbed_pol,costs_pol);

    assert (res == 0);
    assert (dim_vary < dx);

    char ** keys1 = workspace_get_saved_keys(cp->work);
    char ** keys2 = workspace_get_saved_keys2(cp->work);
    size_t * ind_to_serialize = workspace_get_ind_to_serialize(cp->work);

    size_t nprob_run = 0;
    double ** probs = malloc_dd(N);
    size_t * ind_prob_run = calloc_size_t(N);
    
    for (size_t ii = 0; ii < dx; ii++){
        ind_to_serialize[ii] = fi[ii];
    }
    free(fi); fi = NULL;
    ind_to_serialize[dx] = workspace_get_pi_iter(cp->work);
    ind_to_serialize[dx+1] = workspace_get_pi_subiter(cp->work);


    // For the value
    int * absorbed = calloc_int(N);
    double * costs = calloc_double(N*(2*dx+1));
    fi = calloc_size_t(dx);
    res = mca_get_neighbor_costs(dx,N,x,
                                 bound,vf_iteration,ngrid,xgrid,
                                 fi,&dim_vary,
                                 absorbed,costs);
    free(fi); fi = NULL;
    
    for (size_t ii = 0; ii < N; ii++){
        ind_to_serialize[dim_vary] = ii;
        /* printf("serializing: "); iprint_sz(dx+2,ind_to_serialize); */
        size_t_a_to_char(ind_to_serialize,dx+2,keys1[ii]);
        size_t_a_to_char(ind_to_serialize,dx+1,keys2[ii]);
        
        /* printf("key = %s\n",keys1[ii]); */
        /* printf("key = %s\n",keys2[ii]);         */

        size_t nbytes = 0;
        double * out_stored = htable_get_element(htable_iter,keys1[ii],&nbytes);
        if (out_stored != NULL){
            out[ii] = out_stored[0];            
        }
        else if (absorbed_pol[ii] == 1){ // absorbed cost
            res = dp->boundcost(time,x+ii*dx,out+ii);            
            assert (res == 0);
            htable_add_element(htable,keys1[ii],out+ii,1);
            param->niter_evals++; // ran it
            param->niter_node_evals++;
        }
        else if (absorbed[ii] == -1){
            res = dp->obscost(x+ii*dx,out+ii);
            assert (res == 0);

            htable_add_element(htable,keys1[ii],out+ii,1);
            param->niter_evals++; // ran it
            param->niter_node_evals++;
        }
        else{
            param->niter_evals++; // ran it
            param->niter_node_evals++;
            
            probs[ii] = htable_get_element(htable,keys2[ii],&nbytes);
            if (probs[ii] == NULL){
                param->npol_evals++;
                probs[ii] = calloc_double(2*dx+1 + 2);
                ind_prob_run[nprob_run] = ii;
                nprob_run++;
            }
            else{
                out[ii] = bellmanrhs(dx,du,
                                     probs[ii][2*dx+2], // stagecost
                                     NULL,
                                     dp->discount,
                                     probs[ii],
                                     NULL,
                                     probs[ii][2*dx+1]
                                     ,NULL, // dt
                                     costs+ii*(2*dx+1),NULL);

                htable_add_element(htable,keys1[ii],out+ii,1);
            }                
        }
    }


    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
    #endif
    for (size_t ii = 0; ii < nprob_run; ii++){

        struct Memory mem;
        mem.shared = cp;
        mem.private = ind_prob_run[ii];
                

        double * u = workspace_get_u(cp->work,ii);

        double val;
        int res2 = bellman_optimal(du,u,&val,&mem);

                
        double * drift = workspace_get_drift(cp->work, ind_prob_run[ii]);
        double * diff  = workspace_get_diff (cp->work, ind_prob_run[ii]);

        // build transition probabilities at optimal
        res2 = drift_eval(dp->dyn_drift,time,x+ind_prob_run[ii]*dx,u,drift,NULL);
        res2 = diff_eval (dp->dyn_diff,time,x+ind_prob_run[ii]*dx,u,diff,NULL);
        res = transition_assemble(dx,du,cp->dw,mca->h2,mca->t,
                                  drift,NULL,diff,NULL,
                                  probs[ind_prob_run[ii]],
                                  NULL,probs[ind_prob_run[ii]] + 2*dx+1, // dt
                                  NULL,NULL);

        res2 = dp->stagecost(time,x+ind_prob_run[ii]*dx,u,
                             probs[ind_prob_run[ii]] + 2*dx+2, NULL);

        assert (res2 == 0);
        out[ind_prob_run[ii]] = bellmanrhs(dx,du,
                                           probs[ind_prob_run[ii]][2*dx+2], // stagecost
                                           NULL,
                                           dp->discount,
                                           probs[ind_prob_run[ii]],
                                           NULL,
                                           probs[ind_prob_run[ii]][2*dx+1]
                                           ,NULL, // dt
                                           costs+ind_prob_run[ii]*(2*dx+1),NULL);
     

    }
    
    for (size_t ii = 0; ii < nprob_run; ii++){
        htable_add_element(htable,keys2[ind_prob_run[ii]],probs[ind_prob_run[ii]],2*dx+1+2);
        htable_add_element(htable,keys1[ii],out+ind_prob_run[ii],1);
        free(probs[ind_prob_run[ii]]); probs[ind_prob_run[ii]] = NULL;
    }
    free(ind_prob_run); ind_prob_run = NULL;
    free(probs); probs = NULL;
    free(absorbed); absorbed = NULL;
    free(costs); costs = NULL;
    return 0;
}


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


struct C3Control
{
    size_t dx;
    size_t du;
    size_t dw;
    
    size_t * ngrid;
    double ** xgrid;
    double * h;
    double hmin;
    
    struct Boundary * bound;
    struct MCAparam * mca;
    struct DPparam * dp;
    struct Workspace * work;
    
    struct ValueF * policy_sim;
    struct c3Opt * opt_sim;
    void (*transform_sim)(size_t,const double*,double*);
    double * prevpol;
};

struct C3Control *
c3control_create(size_t dx, size_t du, size_t dw,
                 double * lb, double * ub,
                 size_t * ngrid, double discount)
{
    struct C3Control * c3c = malloc(sizeof(struct C3Control));
    c3c->dx = dx;
    c3c->du = du;
    c3c->dw = dw;

    c3c->ngrid = ngrid;
    c3c->xgrid = malloc(dx * sizeof(double *));
    c3c->h     = calloc_double(dx);
    c3c->hmin  = ub[0] - lb[0];
    size_t maxn = ngrid[0];
    for (size_t ii = 0; ii < dx; ii++){
        c3c->xgrid[ii] = linspace(lb[ii],ub[ii],ngrid[ii]);
        c3c->h[ii] = c3c->xgrid[ii][1] - c3c->xgrid[ii][0];
        if (c3c->h[ii] < c3c->hmin){
            c3c->hmin = c3c->h[ii];
        }
        if (ngrid[ii] > maxn){
            maxn = ngrid[ii];
        }
    }
    c3c->bound = boundary_alloc(dx,lb,ub);
    c3c->mca   = mca_param_create(dx,du);
    mca_add_grid_refs(c3c->mca,c3c->ngrid,c3c->xgrid,c3c->hmin,c3c->h);
    c3c->dp    = dp_param_create(dx,du,dw,discount);
    dp_param_add_boundary(c3c->dp, c3c->bound);
    c3c->work = workspace_alloc(dx,du,dw,maxn);

    c3c->policy_sim = NULL;
    c3c->opt_sim = NULL;
    c3c->transform_sim = NULL;
    c3c->prevpol = NULL;
    return c3c;
}

void c3control_destroy(struct C3Control * c3c)
{
    if (c3c != NULL){
        boundary_free(c3c->bound);   c3c->bound = NULL;
        mca_param_destroy(c3c->mca); c3c->mca   = NULL;
        dp_param_destroy(c3c->dp);   c3c->dp    = NULL;
        for (size_t ii = 0; ii < c3c->dx; ii++){
            free(c3c->xgrid[ii]); c3c->xgrid[ii] = NULL;
        }
        workspace_free(c3c->work);  c3c->work = NULL;
        free(c3c->prevpol); c3c->prevpol = NULL;
        free(c3c->xgrid);   c3c->xgrid = NULL;
        free(c3c->h);       c3c->h     = NULL;
        free(c3c);          c3c = NULL;
    }
}

size_t * c3control_get_ngrid(struct C3Control * c3c)
{
    if (c3c == NULL){
        return NULL;
    }
    return c3c->ngrid;
}

double ** c3control_get_xgrid(struct C3Control * c3c)
{
    if (c3c == NULL){
        return NULL;
    }
    return c3c->xgrid;
}

void c3control_add_policy_sim(struct C3Control * c3c, struct ValueF * pol,
                              struct c3Opt * opt_sim,
                              void (*transform)(size_t, const double *, double*))
{
    assert (c3c != NULL);
    c3c->policy_sim = pol;
    c3c->opt_sim = opt_sim;
    c3c->transform_sim = transform;
}

/**********************************************************//**
    Set the external boundary type of a particular dimension
**************************************************************/
void c3control_set_external_boundary(struct C3Control * c3c, size_t dim,
                                     char * type)
{
    assert (c3c != NULL);
    assert (c3c->bound != NULL);
    boundary_external_set_type(c3c->bound,dim,type);
}

void c3control_add_obstacle(struct C3Control * c3c, double * center, double * widths)
{
    assert (c3c != NULL);
    assert (c3c->bound != NULL);
    boundary_add_obstacle(c3c->bound, center, widths);
}

void c3control_add_drift(struct C3Control * c3c, int (*b)(double,const double*,const double*,
                                                          double*,double*,void*),
                         void * args)
{
    assert (c3c != NULL);
    assert (c3c->dp != NULL);
    dp_param_add_drift(c3c->dp,b,args);
}

void c3control_add_diff(struct C3Control * c3c, int (*s)(double,const double*,const double*,
                                                         double*,double*,void*),
                        void * sargs)
{
    assert (c3c != NULL);
    assert (c3c->dp != NULL);
    dp_param_add_diff(c3c->dp,s,sargs);
}

void c3control_add_stagecost(struct C3Control * c3c,
                             int (*stagecost)(double,const double*,
                                              const double*,double*,double*))
{
    assert (c3c != NULL);
    assert (c3c->dp != NULL);
    dp_param_add_stagecost(c3c->dp,stagecost);
}

void c3control_add_boundcost(struct C3Control * c3c,
                             int (*boundcost)(double,const double*,double*))
{
    assert (c3c != NULL);
    assert (c3c->dp != NULL);
    dp_param_add_boundcost(c3c->dp,boundcost);
}

void c3control_add_obscost(struct C3Control * c3c, int (*obscost)(const double*,double*))
{
    assert (c3c != NULL);
    assert (c3c->dp != NULL);
    dp_param_add_obscost(c3c->dp,obscost);
}


int c3control_policy_eval(struct C3Control * c3c, double t,
                          const double * x, double * u)
{
    assert (c3c != NULL);
    assert (c3c->policy_sim != NULL);

    size_t dx = c3c->dx;
    size_t * ngrid = c3c->ngrid;
    double ** xgrid = c3c->xgrid;


    int * absorbed = workspace_get_absorbed(c3c->work,0);
    double * costs = workspace_get_costs(c3c->work,0);

    struct Boundary * bound = c3c->bound;
    struct ValueF * policy_sim = c3c->policy_sim;
    struct c3Opt * opt = c3c->opt_sim;
    assert (policy_sim != NULL);
    assert (opt != NULL);

    struct ControlParams * cp = control_params_create(c3c->dx,c3c->dw,
                                                      c3c->dp,c3c->mca,
                                                      c3c->work,opt);

    double val;
    int res = mca_get_neighbor_node_costs(dx, x, bound, policy_sim,
                                          ngrid, xgrid, absorbed, costs);


    assert (res == 0);
    if (c3c->prevpol == NULL){
        c3c->prevpol = calloc_double(c3c->du);
    }

    for (size_t ii = 0; ii < c3c->du; ii++){
        u[ii] = c3c->prevpol[ii];
    }

    control_params_add_time_and_states(cp,t,1,x);

    struct Memory mem;
    mem.shared = cp;
    mem.private = 0;
    int res2 = bellman_optimal(c3c->du,u,&val,&mem);
    /* printf("val = %G u choose = ",val); dprint(c3c->du,u); */
    assert (res2 == 0);
    /* exit(1); */
    for (size_t ii = 0; ii < c3c->du; ii++){
        c3c->prevpol[ii] = u[ii];
    }

    control_params_destroy(cp); cp = NULL;
    return 0;
}

int c3control_controller(double t,const double * x, double * u, void * args)
{

    struct C3Control * c3c = args;
    if (c3c->transform_sim == NULL){
        return c3control_policy_eval(args,t,x,u);
    }
    else{
        double * xin = calloc_double(c3c->dx);
        c3c->transform_sim(c3c->dx,x,xin);
        int out = c3control_policy_eval(args,t,xin,u);
        free(xin); xin = NULL;
        return out;
    }
    /* return 0; */
}

struct ValueF * c3control_step_vi(struct C3Control * c3c, struct ValueF * vf,
                                  struct ApproxArgs * apargs,
                                  struct c3Opt * opt,
                                  int verbose,
                                  size_t * nevals)
{
    assert (c3c != NULL);
    assert (c3c->ngrid != NULL);
    assert (c3c->xgrid != NULL);
    
    size_t dx = c3c->dx;
    size_t * ngrid = c3c->ngrid;
    double ** xgrid = c3c->xgrid;


    struct ControlParams * cp = control_params_create(c3c->dx,c3c->dw,c3c->dp,c3c->mca,c3c->work,opt);
    struct VIparam * vi = vi_param_create(1e-10);
    vi_param_add_cp(vi,cp);
    vi_param_add_value(vi,vf);
    workspace_increment_vi_iter(c3c->work);
    vi->nstate_evals = 0;
    vi->nnode_evals = 0;
    vi->time_in_loop = 0.0;
    
    struct ValueF * next = valuef_interp(dx,bellman_vi,vi,ngrid,xgrid,vf,apargs,verbose);

    *nevals = vi->nnode_evals;
    if (verbose > 0){
        printf("time in iteration = %G\n",vi->time_in_loop);
    }
    
    vi_param_destroy(vi); vi = NULL;
    control_params_destroy(cp); cp = NULL;
    
    return next;
}

struct ValueF * c3control_step_pi(struct C3Control * c3c, struct ValueF * vf,
                                  struct PIparam * poli,
                                  struct ApproxArgs * apargs,
                                  struct c3Opt * opt,
                                  int verbose,
                                  size_t * niter_evals)
{
    assert (c3c != NULL);
    assert (c3c->ngrid != NULL);
    assert (c3c->xgrid != NULL);
    assert (vf != NULL);
    assert (poli != NULL);
    assert (opt != NULL);
    assert (apargs != NULL);
    assert (c3c->mca != NULL);
    assert (c3c->dp != NULL);
    
    size_t dx = c3c->dx;
    size_t * ngrid = c3c->ngrid;
    double ** xgrid = c3c->xgrid;
 

    struct ControlParams * cp = control_params_create(c3c->dx,c3c->dw,c3c->dp,c3c->mca,c3c->work,opt);

    pi_param_add_cp(poli,cp);
    pi_param_add_value(poli,vf);
    // note this denotes change of cost function and thus this increment
    // changes much faster than pi_subiter.
    // pi_iter is incremented only when a new policy is created
    // in pi_solve
    workspace_increment_pi_subiter(c3c->work); 
    

    /* printf("interp\n"); */
    poli->niter_evals = 0;
    poli->niter_node_evals = 0;

    /* printf("Begin interpolation for PI\n"); */
    /* struct ValueF * next = valuef_interp(dx,bellman_pi,poli,ngrid,xgrid, */
    /*                                      start,apargs,verbose); */
    struct ValueF * next = valuef_interp(dx,bellman_pi,poli,ngrid,xgrid,vf,apargs,verbose);

    /* printf("Completed interpolation for PI\n"); */
    *niter_evals = poli->niter_node_evals;

    control_params_destroy(cp); cp = NULL;

    return next;
}


struct ValueF *
c3control_init_value(struct C3Control * c3c,
                     int (*f)(size_t,const double *,double*,void*),void * args,
                     struct ApproxArgs * aargs, int verbose)
{
    assert (c3c != NULL);
    assert (c3c->ngrid != NULL);
    assert (c3c->xgrid != NULL);
    
    struct ValueF * vf = valuef_interp(c3c->dx,f,args,c3c->ngrid,c3c->xgrid,
                                       NULL, aargs,verbose);


    return vf;
        
}

struct ValueF * c3control_vi_solve(struct C3Control * c3c,
                                   size_t maxiter, double abs_conv_tol,
                                   struct ValueF * vo,
                                   struct ApproxArgs * apargs,
                                   struct c3Opt * opt,
                                   int verbose, 
                                   struct Diag ** diag)
{

    struct ValueF * start = valuef_copy(vo);
    workspace_reset_vi_htable(c3c->work);
    for (size_t ii = 0; ii < maxiter; ii++){

        // memory
        if (ii % 1000 == 0){ // precaution against too much growth in memory usage
            workspace_reset_vi_htable(c3c->work);
        }
        
        /* printf("ii = %zu\n",ii); */
        assert (start != NULL);
        size_t niter_evals = 0;
        struct ValueF * next = c3control_step_vi(c3c,start,apargs,opt,
                                                 verbose-1,
                                                 &niter_evals);
        /* printf("stepped\n"); */
        double diff = valuef_norm2diff(start,next);
        double norm = valuef_norm(next);
        size_t stot = 1;
        for (size_t jj = 0; jj < c3c->dx; jj++){
            stot *= c3c->ngrid[jj];
        }
        double frac = (double) niter_evals / (double) stot;
        if (verbose > 0){
            if ((ii+1) % 1 == 0){
                printf("\t Value Iteration (%zu\\%zu):\n",ii+1,maxiter);
                printf("\t \t L2 Difference between iterates    = %3.5E\n ",diff);
                printf("\t \t L2 Norm of current value function = %3.5E\n", norm);
                printf("\t \t Relative L2 Cauchy difference     = %3.5E\n", diff/norm);
                printf("\t \t Fraction of states evaluated      = %3.5E\n", frac);
                /* printf("\t \t Ratio: L2 Norm Prev / L2 Norm Cur = %G\n",rat); */
            }
        }
        
        if (diag != NULL){
            size_t * ranks = valuef_get_ranks(next);
            diag_append(diag,ii,1,norm,diff,c3c->dx,ranks,frac);
        }

        
        valuef_destroy(start); start = NULL;
        start = valuef_copy(next);
        valuef_destroy(next); next = NULL;
        norm = valuef_norm(start); 
        if (diff < abs_conv_tol){
            break;
        }
    }
    return start;
}


struct ValueF * c3control_pi_solve(struct C3Control * c3c,
                                   size_t maxiter, double abs_conv_tol,
                                   struct ValueF * policy,
                                   struct ApproxArgs * apargs,
                                   struct c3Opt * opt,
                                   int verbose, struct Diag ** diag)
{

    struct ValueF * start = valuef_copy(policy);
    struct PIparam * poli = pi_param_create(1e-10,policy);
    workspace_increment_pi_iter(c3c->work);
    workspace_reset_pi_prob_htable(c3c->work);
    workspace_reset_pi_htable(c3c->work);
    
    /* double dd = valuef_norm2diff(start,policy); */
    /* printf("diff = %G\n",dd); */
    for (size_t ii = 0; ii < maxiter; ii++){
        /* printf("ii = %zu\n",ii); */
        assert (start != NULL);
        size_t niter_evals = 0;
        struct ValueF * next = c3control_step_pi(c3c,start,poli,apargs,opt,
                                                 verbose-1,
                                                 &niter_evals);
        /* printf("stepped\n"); */
        double diff = valuef_norm2diff(start,next);
        double norm = valuef_norm(next);
        size_t stot = 1;
        for (size_t jj = 0; jj < c3c->dx; jj++){
            stot *= c3c->ngrid[jj];
        }
        double frac = (double) niter_evals / (double) stot;

        if (verbose > 0){
            if ((ii+1) % 1 == 0){
                printf("\t POLICY ITERATION (%zu\\%zu):\n",ii+1,maxiter);
                printf("\t \t L2 Difference between iterates    = %3.5E\n ",diff);
                printf("\t \t L2 Norm of current value function = %3.5E\n", norm);
                printf("\t \t Relative L2 Cauchy difference     = %3.5E\n", diff/norm);
                printf("\t \t Fraction of states evaluated      = %3.5E\n", frac);
                /* printf("\t \t Ratio: L2 Norm Prev / L2 Norm Cur = %G\n",rat); */
            }
        }
        
        if (diag != NULL){
            /* printf("frac = %G\n",frac); */
            size_t * ranks = valuef_get_ranks(next);
            diag_append(diag,ii,0,norm,diff,c3c->dx,ranks,frac);
            /* printf("appended\n"); */
        }
        
        valuef_destroy(start); start = NULL;
        start = valuef_copy(next);
        valuef_destroy(next); next = NULL;
        norm = valuef_norm(start); 
        /* printf("start norm = %G\n",norm); */

        
        if (diff < abs_conv_tol){
            break;
        }

    }
    pi_param_destroy(poli); poli = NULL;
    return start;
}

struct Diag
{
    size_t iter;
    int type; // 0 for pi, 1 for vi;
    double norm;
    double abs_diff;
    size_t dim;
    size_t * ranks;
    double frac;

    struct Diag * next;
};

void diag_destroy(struct Diag ** head)
{
    if (*head != NULL){
        struct Diag * current = *head;
        struct Diag * next;
        while (current != NULL){
            next = current->next; current->next = NULL;
            free(current->ranks); current->ranks = NULL;
            free(current); current = NULL;
            current = next;
        }
        *head = NULL;
    }
}

struct Diag * diag_create(size_t iter, int type, double norm,
                          double abs_diff, size_t dim, size_t * ranks,
                          double frac)
{
    struct Diag * diag = malloc(sizeof(struct Diag));
    diag->iter = iter;
    diag->type = type;
    diag->norm = norm;
    diag->abs_diff = abs_diff;
    diag->dim = dim;
    diag->ranks = calloc_size_t(dim+1);
    for (size_t ii = 0; ii < dim; ii++){
        diag->ranks[ii] = ranks[ii];
    }
    diag->ranks[dim] = ranks[dim];
    diag->frac= frac;
    diag->next = NULL;
    return diag;
}

void diag_append(struct Diag ** diag,
                 size_t iter, int type, double norm,
                 double abs_diff, size_t dim, size_t * ranks,
                 double frac)
{
    struct Diag * current = *diag;
    if (current == NULL){
        /* printf("create\n"); */
        *diag = diag_create(iter,type,norm,abs_diff,dim,ranks,frac);
        /* printf("created\n"); */
    }
    else{
        while (current->next != NULL){
            current = current->next;
        }
        current->next = diag_create(iter,type,norm,abs_diff,dim,ranks,frac);
    }
}

void diag_print(struct Diag * head, FILE * fp)
{
    struct Diag * current = head;
    while ( current != NULL){
        double avg_rank = 0.0;
        size_t maxrank = 0;
        for (size_t ii = 1; ii < current->dim; ii++){
            if (current->ranks[ii] > maxrank){
                maxrank = current->ranks[ii];
            }
            avg_rank += (double)(current->ranks[ii]);
        }
        /* printf("ranks = "); iprint_sz(current->dim+1,current->ranks); */
        avg_rank /= (double)(current->dim-1);
        fprintf(fp,"%zu %d %3.15G %3.15G %3.15G %zu %3.15G \n",
                current->iter, 
                current->type, 
                current->norm,
                current->abs_diff,
                avg_rank,
                maxrank,
                current->frac);
        current = current->next;
    }
}

int diag_save(struct Diag * head, char * filename)
{
    FILE * fp =  fopen(filename, "w");
    if (fp == NULL){
        fprintf(stderr, "cat: can't open %s\n", filename);
        return 1;
    }

    diag_print(head,fp);
    
    fclose(fp);
    return 0;
}
