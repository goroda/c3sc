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
        /* printf("dtgrad = "); dprint(du,dtgrad); */
        /* printf("stagegrad = "); dprint(du,stage_grad); */
        /* printf("stagecost = %G\n",stage_cost); */
        /* printf("dt = %G\n",dt); */
        /* printf("ebt = %G\n",ebt); */
        for (size_t jj = 0; jj < du; jj++){
            /* printf("jj = %zu\n",jj); */
            // derivative of stage cost
            grad[jj] = stage_grad[jj] * dt + dtgrad[jj]*stage_cost ;
            grad[jj] += (-discount) * dtgrad[jj] * ebt * cost_to_go;
            /* printf("grad[jj] til noow = %G\n",grad[jj]); */
            for (size_t ii = 0; ii < 2*dx+1; ii++){
                /* printf("ii = %zu, cost=%G\n",ii,cost[ii]); */
                /* printf("prob grad = %G\n",prob_grad[jj*(2*dx+1)+ii]);  */
                grad[jj] += ebt * prob_grad[ii*du + jj] * cost[ii];
                /* printf("grad[jj] is now %G\n",grad[jj]); */

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

    double dt;
    double * grad_dt; // (du,1)

    double * prob; // (2dx+1,1)
    double * grad_prob; //(du * (2dx+1),1)
    
    double * workspace; //(du,1);
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
    mca->grad_dt = calloc_double(du);
    mca->prob = calloc_double(2*dx+1);
    mca->grad_prob = calloc_double(du * (2*dx+1));
    mca->workspace = calloc_double(du);

    mca->ngrid = NULL;
    mca->xgrid = NULL;
    mca->hvec = NULL;
    mca->hmin = 0;
    
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
}

/**********************************************************//**
    Destroy memory allocated to MCA structure
                                                           
    \param[in,out] mca   - MCA structure to modify    
**************************************************************/
void mca_param_destroy(struct MCAparam * mca)
{
    if (mca != NULL){
        free(mca->grad_dt); mca->grad_dt = NULL;
        free(mca->prob); mca->prob = NULL;
        free(mca->grad_prob); mca->grad_prob = NULL;
        free(mca->workspace); mca->workspace = NULL;
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
    double * drift;
    double * grad_drift;
    double * diff;
    double * grad_diff;

    struct Boundary * bound;
    
    // cost functions
    double discount; // discount factor
    int (*stagecost)(double,const double*,const double*,double*,double*);
    double * grad_stage;
    int (*boundcost)(double,const double*,double*);
    int (*obscost)(const double*,double*);
};

struct DPparam * dp_param_create(size_t dx, size_t du, size_t dw, double discount)
{
    struct DPparam * dp = malloc(sizeof(struct DPparam));
    assert (dp != NULL);
    
    dp->dyn_drift = drift_alloc(dx,du);
    dp->dyn_diff = diff_alloc(dx,du,dw);

    dp->drift = calloc_double(dx);
    dp->grad_drift = calloc_double(dx*du);
    dp->diff = calloc_double(dx*dw);
    /* printf("allocated diff\n"); */
    /* dprint2d_col(dx,dw,dp->diff); */
    dp->grad_diff = calloc_double(dx*dw*du);

    dp->bound = NULL;
    
    dp->discount = discount;
    dp->grad_stage = calloc_double(du);

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
        free(dp->drift); dp->drift = NULL;
        free(dp->grad_drift); dp->grad_drift = NULL;
        free(dp->diff); dp->diff = NULL;
        free(dp->grad_diff); dp->grad_diff = NULL;
        free(dp->grad_stage); dp->grad_stage = NULL;
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
    const double * x;
    int absorbed;
    const double * costs;
    struct DPparam * dp;
    struct MCAparam * mca;
    struct c3Opt * opt;

    int res_last_grad; // result of the last gradient evaluation
    
};

struct ControlParams * control_params_create(size_t dx, size_t dw, struct DPparam * dp,
                                             struct MCAparam * mca, struct c3Opt * opt)
{
    struct ControlParams * c = malloc(sizeof(struct ControlParams));
    assert (c != NULL);
    c->dx = dx;
    c->dw = dw;
    c->dp = dp;
    c->mca = mca;
    c->opt = opt;
    c->res_last_grad = 0;
    return c;
}

void control_params_add_state_info(struct ControlParams * cp,
                                   double time, const double * x, int absorbed,
                                   const double * costs)
{
    assert (cp != NULL);
    cp->time = time;
    cp->x = x;
    cp->absorbed = absorbed;
    cp->costs = costs;
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
double bellman_control(size_t du, double * u, double * grad_u, void * args)
{
    struct ControlParams * param = args;

    // unpack
    struct DPparam * dp = param->dp;
    struct MCAparam * mca = param->mca;
    const double * costs = param->costs;
    const double * x = param->x;
    double time = param->time;
    size_t dx = param->dx;
    size_t dw = param->dw;
    int absorbed = param->absorbed;
        
    //get drift and diffusion
    int res;
    if (grad_u != NULL){
        for (size_t ii = 0; ii < du; ii++){
            grad_u[ii] = 0.0;
        }
        res = drift_eval(dp->dyn_drift,time,x,u,
                         dp->drift,dp->grad_drift);
        /* printf("drift = "); dprint(dx,dp->drift); */
        assert (res == 0);
        res = diff_eval(dp->dyn_diff,time,x,u,
                        dp->diff,dp->grad_diff);
    }
    else{
        /* printf("x = ");dprint(dx,x); */
        res = drift_eval(dp->dyn_drift,time,x,u,
                         dp->drift,NULL);
        assert (res == 0);
        /* printf("drift = "); dprint(dx,dp->drift); */
        res = diff_eval(dp->dyn_diff,time,x,u,
                        dp->diff,NULL);
        /* printf("diffusion = \n"); */
        /* dprint2d_col(dx,dw,dp->diff); */
        /* exit(1); */
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

            param->res_last_grad = res;
            /* assert (res == 0); */
            val = bellmanrhs(dx,du,stage_cost,dp->grad_stage,dp->discount,
                             mca->prob,mca->grad_prob,mca->dt,mca->grad_dt,
                             costs, grad_u);
        }
        else{
            /* printf("here\n"); */
            res = dp->stagecost(time,x,u,&stage_cost,NULL);
            /* printf("stagecost = %G\n",stage_cost); */
            /* printf("drift = "); dprint(dx,dp->drift); */
            /* printf("diff =  "); dprint2d_col(dx,dw,dp->diff); */
            assert (res == 0);
            /* printf("assemble transition\n"); */
            res = transition_assemble(dx,du,dw,mca->hmin,mca->hvec,
                                      dp->drift,NULL,dp->diff,NULL,mca->prob,
                                      NULL,&(mca->dt),NULL,NULL);
            /* printf("bellmanrhs \n"); */
            /* printf("mca->prob = \n"); */
            /* dprint(2*dx+1,mca->prob); */
            assert (res == 0);
            val = bellmanrhs(dx,du,stage_cost,NULL,dp->discount,mca->prob,NULL,mca->dt,NULL,
                             costs,NULL);
        }
    }
    else if (absorbed == 1){ // absorbed cost
        /* printf("asborbed\n"); */
        res = dp->boundcost(time,x,&val);
        /* printf("val = %G\n",val); */
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

/**********************************************************//**
    Find the optimal control

    \param[in]     du   - number of control variables
    \param[in,out] u    - control
    \param[in]     val  - value of bellman function at optimal control
    \param[in]     args - pointer to problem parameters (ControlParams)

    \return 0 if successful, otherwise not
**************************************************************/
int bellman_optimal(size_t du, double * u, double * val, void * arg)
{
    assert (arg != NULL);
    struct ControlParams * param = arg;

    // unpack
    struct c3Opt * opt = param->opt;
    assert (opt != NULL);
    double * lbu = c3opt_get_lb(opt);
    double * ubu = c3opt_get_ub(opt);
    double * umin = calloc_double(du);
    double * ucurr = calloc_double(du);
    double minval = 0.0;
    size_t nrand = 10; // note this is changed if du = 1 or du = 3;
    size_t npert = 2;
    double valtemp;
    int justrand = 1;

    c3opt_add_objective(opt,&bellman_control,param);
    // first do zero;
    valtemp = bellman_control(du,ucurr,NULL,arg);
    minval = valtemp; 
    memmove(umin,ucurr,du*sizeof(double));


    if (justrand == 1){
        // now compare with random samples
        nrand = 200;
        for (size_t jj = 0; jj < nrand; jj++){
            for (size_t kk = 0; kk < du; kk++){
                ucurr[kk] = randu()*(ubu[kk]-lbu[kk]) + lbu[kk];
            }
            valtemp = bellman_control(du,ucurr,NULL,arg);
            if (valtemp < minval){
                minval = valtemp;
                memmove(umin,ucurr,du*sizeof(double));
            }
        }
    }
    else{
        // then start at each corner
        /* printf("corners\n"); */
        if (du == 1){
            nrand = 0;
            double distx = (ubu[0]-lbu[0])/4.0;
            ucurr[0] = lbu[0] + distx;
            c3opt_minimize(opt,ucurr,&valtemp);
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
                /* printf("\t updating current min!\n"); */
                minval = valtemp;
                memmove(umin,ucurr,du*sizeof(double));
            }
        }
        else if (du == 3){
            nrand = 0;
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
                        if (valtemp < minval){
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
                            if (valtemp < minval){
                                /* printf("\t updating current min!\n"); */
                                minval = valtemp;
                                memmove(umin,ut,du*sizeof(double));
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

    double convergence;
};

struct VIparam * vi_param_create(double convergence)
{
    struct VIparam * vi = malloc(sizeof(struct VIparam));
    assert (vi != NULL);
    vi->convergence = convergence;

    vi->cp = NULL;
    vi->vf = NULL;

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
}

int bellman_vi(size_t N, const double * x, double * out, void * arg)
{
    struct VIparam * param = arg;
    struct ControlParams * cp = param->cp;
    struct MCAparam * mca = cp->mca;
    struct DPparam * dp = cp->dp;
    size_t dx = mca->dx;
    size_t du = mca->du;
    size_t * ngrid = mca->ngrid;
    double ** xgrid = mca->xgrid;
    
    struct ValueF * vf = param->vf;
    struct Boundary * bound = dp->bound;

    int * absorbed = calloc_int(N);
    double * costs = calloc_double(N*(2*dx+1));
    double * u = calloc_double(du);

    /* printf("lets go!\n"); */
    /* for (size_t ii = 0; ii < N; ii++){ */
    /*     dprint(dx,x+ii*dx); */
    /* } */
    /* printf("get neighbor cost\n"); */
    int res = mca_get_neighbor_costs(dx,N,x,
                                     bound,vf,ngrid,xgrid,
                                     absorbed,costs);
    /* printf("got it\n"); */
    assert (res == 0);

    double time = 0.0;

    for (size_t ii = 0; ii < N; ii++){
        control_params_add_state_info(cp,time,x+ii*dx,absorbed[ii],costs+ii*(2*dx+1));
        int res2 = bellman_optimal(du,u,out+ii,cp);
        assert (res2 == 0);
    }

    free(absorbed); absorbed = NULL;
    free(costs); costs = NULL;
    free(u); u = NULL;

    return 0;
    
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
