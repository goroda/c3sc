#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "c3.h"
#include "util.h"
#include "nodeutil.h"
#include "dynamics.h"
#include "hashgrid.h"

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

struct ControlParams * control_params_create(size_t dx, size_t dw,
                                             struct DPparam * dp,
                                             struct MCAparam * mca, 
                                             struct c3Opt * opt)
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
            /* Printf("Stagecost = %G\n",stage_cost); */
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

    \param[in]     du  - number of control variables
    \param[in,out] u   - control
    \param[in]     val - value of bellman function at optimal control
    \param[in]     arg - pointer to problem parameters (ControlParams)

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
    int justrand = 0;

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
        c3opt_minimize(opt,umin,val);
        
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
            nrand = 10;
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
    struct HTable * htable;

    double convergence;
};

struct VIparam * vi_param_create(double convergence)
{
    struct VIparam * vi = malloc(sizeof(struct VIparam));
    assert (vi != NULL);
    vi->convergence = convergence;
    vi->htable = htable_create(1000000);
    vi->cp = NULL;
    vi->vf = NULL;

    return vi;
}

void vi_param_destroy(struct VIparam * vi)
{
    if (vi != NULL){
        htable_destroy(vi->htable); vi->htable = NULL;
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

    // create a new htable for the new value function
    htable_destroy(vi->htable); vi->htable = NULL;
    vi->htable = htable_create(1000000);
}

int bellman_vi(size_t N, const double * x, double * out, void * arg)
{
    struct VIparam * param = arg;
    struct ControlParams * cp = param->cp;
    assert (cp != NULL);
    struct MCAparam * mca = cp->mca;
    struct DPparam * dp = cp->dp;
    size_t dx = mca->dx;
    size_t du = mca->du;
    size_t * ngrid = mca->ngrid;
    double ** xgrid = mca->xgrid;
    
    struct ValueF * vf = param->vf;
    struct HTable * htable = param->htable;
    assert (vf != NULL);
    struct Boundary * bound = dp->bound;

    int * absorbed = calloc_int(N);
    double * costs = calloc_double(N*(2*dx+1));
    double * u = calloc_double(du);

    /* printf("lets go!\n"); */
    /* for (size_t ii = 0; ii < N; ii++){ */
    /*     dprint(dx,x+ii*dx); */
    /* } */
    /* printf("get neighbor cost\n"); */
    size_t * fi = calloc_size_t(dx);
    size_t dim_vary;
    int res = mca_get_neighbor_costs(dx,N,x,
                                     bound,vf,ngrid,xgrid,
                                     fi,&dim_vary,
                                     absorbed,costs);
    assert (res == 0);
    size_t * ind_to_serialize = calloc_size_t(dx+1);
    for (size_t ii = 0; ii < dx; ii++){
        if (ii != dim_vary){
            ind_to_serialize[ii] = fi[ii];
        }
    }
    ind_to_serialize[dx] = dim_vary;
    free(fi); fi = NULL;
    char * key = size_t_a_to_char(ind_to_serialize,dx+1);
    free(ind_to_serialize); ind_to_serialize = NULL;
    /* free(fi); fi = NULL; */

    size_t nbytes = 0;
    double * out_stored = htable_get_element(htable,key,&nbytes);

    if (out_stored == NULL){
        /* printf("got it\n"); */
        double time = 0.0;
        for (size_t ii = 0; ii < N; ii++){
            control_params_add_state_info(cp,time,x+ii*dx,absorbed[ii],
                                          costs+ii*(2*dx+1));
            int res2 = bellman_optimal(du,u,out+ii,cp);
            assert (res2 == 0);
        }
        htable_add_element(htable,key,out,N * sizeof(double));
    }
    else{
        for (size_t ii = 0; ii < N; ii++){
            out[ii] = out_stored[ii];
        }
        free(key); key = NULL;
    }

    free(absorbed); absorbed = NULL;
    free(costs); costs = NULL;
    free(u); u = NULL;

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
    struct HTable * htable_iter; // store evaluations for fixed iteration

    struct ValueF * vf_policy;
    struct HTable * htable; // store probability transitions, stage costs, etc associated with policy
    
    double convergence;
};

struct PIparam * pi_param_create(double convergence, struct ValueF * policy)
{
    struct PIparam * poli = malloc(sizeof(struct PIparam));
    assert (poli != NULL);
    poli->convergence = convergence;
    poli->vf_policy = policy;
    poli->htable = htable_create(100000);
        
    poli->cp = NULL;
    poli->vf_iteration = NULL;
    poli->htable_iter = NULL;

    return poli;
}

void pi_param_destroy(struct PIparam * poli)
{
    if (poli != NULL){
        htable_destroy(poli->htable); poli->htable = NULL;
        htable_destroy(poli->htable_iter); poli->htable_iter = NULL;
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

    htable_destroy(poli->htable_iter); 
    poli->htable_iter = NULL;
    poli->htable_iter = htable_create(1000000);
}

int bellman_pi(size_t N, const double * x, double * out, void * arg)
{
    struct PIparam * param = arg;
    assert (param->cp != NULL);
    struct ControlParams * cp = param->cp;
    struct MCAparam * mca = cp->mca;
    struct DPparam * dp = cp->dp;
    size_t dx = mca->dx;
    size_t du = mca->du;
    size_t * ngrid = mca->ngrid;
    double ** xgrid = mca->xgrid;
    
    struct ValueF * vf_policy    = param->vf_policy;
    struct ValueF * vf_iteration = param->vf_iteration;
    struct Boundary * bound = dp->bound;

    struct HTable * htable = param->htable;
    struct HTable * htable_iter = param->htable_iter;
    
    // handle the policy
    double time = 0.0;


    int * absorbed_pol = calloc_int(N);
    double * costs_pol = calloc_double(N*(2*dx+1));

    size_t * fi = calloc_size_t(dx);
    size_t dim_vary;
    int res = mca_get_neighbor_costs(dx,N,x,
                                     bound,vf_policy,ngrid,xgrid,
                                     fi,&dim_vary,
                                     absorbed_pol,costs_pol);

    size_t * ind_to_serialize = calloc_size_t(dx+1);
    for (size_t ii = 0; ii < dx; ii++){
        if (ii != dim_vary){
            ind_to_serialize[ii] = fi[ii];
        }
    }
    ind_to_serialize[dx] = dim_vary;
    free(fi); fi = NULL;
    char * key1 = size_t_a_to_char(ind_to_serialize,dx+1);
    char * key2 = size_t_a_to_char(ind_to_serialize,dx+1);
    char * key3 = size_t_a_to_char(ind_to_serialize,dx+1);
    free(ind_to_serialize); ind_to_serialize = NULL;
    


    // see if already computed the output for this function
    size_t nbytes1 = 0;
    double * out_stored = htable_get_element(htable_iter,key1,&nbytes1);
    /* double * out_stored = NULL; */
    if (out_stored == NULL){
        size_t nbytes = 0;
        // see if already computed transition probabilities for this fiber
        double * probs = htable_get_element(htable,key1,&nbytes);
        int probs_alloc = 0;
        /* probs = NULL; */
        if (probs == NULL){
            /* printf("probs is null key=%s\n",key); */
            size_t nprobs = N*(2*dx+1) + 2*N; // last N are stage costs and N before that are dts
            probs = calloc_double(nprobs); // last N are dts
            probs_alloc = 1;

            double * u = calloc_double(du);
            for (size_t ii = 0; ii < N; ii++){
                control_params_add_state_info(cp,time,x+ii*dx,
                                              absorbed_pol[ii],
                                              costs_pol+ii*(2*dx+1));
                double val;
                int res2 = bellman_optimal(du,u,&val,cp);

                // build transition probabilities
                assert (res2 == 0);
                res2 = drift_eval(dp->dyn_drift,time,x+ii*dx,u,dp->drift,NULL);
                assert (res2 == 0);
                res2 = diff_eval(dp->dyn_diff,time,x+ii*dx,u,dp->diff,NULL);
                res = transition_assemble(dx,du,cp->dw,mca->hmin,mca->hvec,
                                          dp->drift,NULL,dp->diff,NULL,
                                          probs + ii*(2*dx+1),
                                          NULL,probs + N*(2*dx+1)+ii, // dt
                                          NULL,NULL);

                res2 = dp->stagecost(time,x+ii*dx,u,probs + N*(2*dx+1)+N+ii,NULL);
                assert (res2 == 0);
            }
            htable_add_element(htable,key2,probs,nprobs * sizeof(double));
            free(u); u = NULL;
        }
        else{
            free(key2); key2 = NULL;
        }

        // need to compute the control somewhere
        int * absorbed = calloc_int(N);
        double * costs = calloc_double(N*(2*dx+1));
        fi = calloc_size_t(dx);
        res = mca_get_neighbor_costs(dx,N,x,
                                     bound,vf_iteration,ngrid,xgrid,
                                     fi,&dim_vary,
                                     absorbed,costs);

        assert (res == 0);
        free(fi); fi = NULL;

        for (size_t ii = 0; ii < N; ii++){
            // get the optimal control
            // iterate with the optimal control
            if (absorbed[ii] == 0){
                out[ii] = bellmanrhs(dx,du,
                                     *(probs + N*(2*dx+1)+N+ii), // stagecost
                                     NULL,
                                     dp->discount, 
                                     probs+ii*(2*dx+1),NULL,
                                     *(probs+N*(2*dx+1)+ii),NULL, // dt
                                     costs+ii*(2*dx+1),NULL);
            }
            else if (absorbed[ii] == 1){ // absorbed cost
                res = dp->boundcost(time,x+ii*dx,out+ii);
                assert (res == 0);
            }
            else if (absorbed[ii] == -1){
                res = dp->obscost(x+ii*dx,out+ii);
            }
            else{
                fprintf(stderr, "Unrecognized aborbed condition %d\n",absorbed[ii]);
                exit(1);
            }
        }

        if (probs_alloc == 1){
            free(probs); probs = NULL;
        }
        free(absorbed); absorbed = NULL;
        free(costs); costs = NULL;
        
        htable_add_element(htable_iter,key3,out,N * sizeof(double));
    }
    else{
        for (size_t ii = 0; ii < N; ii++){
            out[ii] = out_stored[ii];
        }
        free(key3); key3 = NULL;
        free(key2); key2 = NULL;
    }

    free(key1); key1 = NULL;
    /* free(key2); key2 = NULL; */

    free(absorbed_pol); absorbed_pol = NULL;
    free(costs_pol); costs_pol = NULL;

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
    for (size_t ii = 0; ii < dx; ii++){
        c3c->xgrid[ii] = linspace(lb[ii],ub[ii],ngrid[ii]);
        c3c->h[ii] = c3c->xgrid[ii][1] - c3c->xgrid[ii][0];
        if (c3c->h[ii] < c3c->hmin){
            c3c->hmin = c3c->h[ii];
        }
    }
    c3c->bound = boundary_alloc(dx,lb,ub);
    c3c->mca   = mca_param_create(dx,du);
    mca_add_grid_refs(c3c->mca,c3c->ngrid,c3c->xgrid,c3c->hmin,c3c->h);
    c3c->dp    = dp_param_create(dx,du,dw,discount);
    dp_param_add_boundary(c3c->dp, c3c->bound);

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
        free(c3c->xgrid); c3c->xgrid = NULL;
        free(c3c->h);     c3c->h     = NULL;
        free(c3c); c3c = NULL;
    }
}

void c3control_add_obstacle(struct C3Control * c3c, double * center, double * widths)
{
    assert (c3c != NULL);
    assert (c3c->bound == NULL);
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
