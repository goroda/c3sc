#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "c3.h"

#include "util.h"
#include "simulate.h"
#include "control.h"
#include "cost.h"
#include "dp.h"
#include "tensmarkov.h"

/** \struct C3SC
 *  \brief Stochastic control problem
 *  \var C3SC::type
 *  problem type
 *  \var C3SC::dx
 *  dimension of state
 *  \var C3SC::du
 *  dimension of control
 *  \var C3SC::dw
 *  dimension of noise
 *  \var C3SC::dpih
 *  infinite horizon problem
 *  \var C3SC::dpfh
 *  finite horizon problem
 */
struct C3SC
{
    enum SCTYPE type;
    size_t dx;
    size_t du;
    size_t dw;
    double * lbx;
    double * ubx;

    struct DPih * dpih;
    struct DPfh * dpfh;

    // auxilary stuff
    struct Boundary * bound;
    struct Dyn * dyn;
    struct MCA * mca;
    struct Cost * cost;
    struct Policy * pol;
    struct c3Opt * opt;
};

/**********************************************************//**
    Create a stochastic control problem

    \param[in] type  - IH (infinite horizon) or FH (finite horizon)
    \param[in] dx    - dimension of state space
    \param[in] du    - dimension of control space
    \param[in] dw    - dimension of noise

    \return stochastic control problem
**************************************************************/
struct C3SC * c3sc_create(enum SCTYPE type, size_t dx, size_t du, size_t dw)
{
    struct C3SC * sc = malloc(sizeof(struct C3SC));
    if (sc == NULL){
        fprintf(stderr,"Cannot create c3sc\n");
        exit(1);
    }
    
    sc->type = type;
    sc->dx = dx;
    sc->du = du;
    sc->dw = dw;
    sc->lbx = NULL;
    sc->ubx = NULL;
    
    sc->bound = NULL;
    sc->dyn = NULL;
    sc->mca = NULL;
    sc->cost = NULL;
    sc->pol = NULL;
    sc->opt = NULL;
    return sc;
}
/**********************************************************//**
    Free memory allocated to stochastic control problem
**************************************************************/
void c3sc_destroy(struct C3SC * sc)
{
    if (sc != NULL){
        assert (sc->type != FH);
        free(sc->lbx); sc->lbx = NULL;
        free(sc->ubx); sc->ubx = NULL;
        dpih_free_deep(sc->dpih); sc->dpih = NULL;
        free(sc); sc = NULL;
    }
}

/**********************************************************//**
    Add bounds on state
**************************************************************/
void
c3sc_set_state_bounds(struct C3SC * sc, double * lb, double * ub)
{
    assert (sc != NULL);
    if (lb != NULL){
        sc->lbx = calloc_double(sc->dx);
        memmove(sc->lbx,lb,sc->dx * sizeof(double));
    }
    if (ub != NULL){
        sc->ubx = calloc_double(sc->dx);
        memmove(sc->ubx,ub,sc->dx * sizeof(double));
    }
}

/**********************************************************//**
    Add dynamics of stochastic cotnrol

    \param[in,out] sc    - stochastic control problem
    \param[in]     b     - drift dynamics
    \param[in]     bargs - drift dynamics arguments
    \param[in]     s     - diffusion dynamics
    \param[in]     sargs - diffusion dynamics arguments
**************************************************************/
void c3sc_add_dynamics(struct C3SC * sc,  
                       int (*b)(double,double*,double*,
                                double*,double*,void*),
                       void * bargs,
                       int (*s)(double,double*,double*,
                                double*,double*,void*),
                       void * sargs)
{
    assert (sc != NULL);
    struct Drift * drift = drift_alloc(sc->dx,sc->du);
    drift_add_func(drift,b,bargs);
    struct Diff * diff = diff_alloc(sc->dx,sc->du,sc->dw);
    diff_add_func(diff,s,sargs);
    sc->dyn = dyn_alloc(drift,diff);
}

/**********************************************************//**
    Add boundary to sc problem

    \param[in,out] sc    - stochastic control problem
    \param[in]     f     - boundary evaluation function
    \param[in]     farg  - additional boundary arguments
**************************************************************/
void c3sc_add_boundary(struct C3SC * sc,
                       int (*f)(double,double *,void *,int*), 
                       void *farg)
{
    assert (sc != NULL);
    sc->bound = boundary_alloc(sc->dx,f,farg);
}

/**********************************************************//**
    Initialize the markov chain approximation method 
    (cost and mca structure)

    \param[in,out] sc - stochastic control problem
    \param[in]     N  - number of nodes in each dimension
**************************************************************/
void c3sc_init_mca(struct C3SC * sc, size_t * N)
{
    assert (sc != NULL);
    if (sc->dyn == NULL){
        fprintf(stderr,"Must add dynamics to c3sc before initializing\n");
        fprintf(stderr,"the Markov Chain approximation algorithm\n");
        exit(1);
    }
    if (sc->bound == NULL){
        fprintf(stderr,"Must add boundary cond. to c3sc before initializing\n");
        fprintf(stderr,"the Markov Chain approximation algorithm\n");
        exit(1);
    }


    sc->cost = cost_alloc(sc->dx,sc->lbx,sc->ubx);
    double ** x = malloc_dd(sc->dx);
    double * h = calloc_double(sc->dx);
    for (size_t ii = 0; ii < sc->dx; ii++){
        x[ii] = linspace(sc->lbx[ii],sc->ubx[ii],N[ii]);
        h[ii] = x[ii][1] - x[ii][0];
    }
    cost_init_discrete(sc->cost,N,x);

    sc->mca = mca_alloc(sc->dx,sc->du,sc->dw,h);
    mca_attach_dyn(sc->mca,sc->dyn);
    mca_attach_bound(sc->mca,sc->bound);
    
    free(h); h = NULL; 
    free_dd(sc->dx,x); x = NULL;
}

/**********************************************************//**
    Attach (by reference )the optimization options and 
    algorithms to the problem

    \param[in,out] sc  - stochastic control problem
    \param[in]     opt - optimization structure
**************************************************************/
void c3sc_attach_opt(struct C3SC * sc, struct c3Opt * opt)
{
    assert(sc != NULL);
    assert(opt != NULL);
    sc->opt = opt;
}


/**********************************************************//**
    Attach (by reference )the optimization options and 
    algorithms to the problem

    \param[in,out] sc   - stochastic control problem
    \param[in]     beta - discount factor
    \param[in]     s    - stage cost
    \param[in]     b    - boundary cost
**************************************************************/
void c3sc_init_dp(struct C3SC * sc, double beta,
                  int (*s)(double,double*,double*,double*,double*),
                  int (*b)(double,double*,double*))
{
    assert(sc != NULL);

    if (sc->dyn == NULL){
        fprintf(stderr,"Must add dynamics to c3sc before initializing\n");
        fprintf(stderr,"the dynamic program\n");
        exit(1);
    }
    if (sc->bound == NULL){
        fprintf(stderr,"Must add boundary cond. to c3sc before initializing\n");
        fprintf(stderr,"the dynamic program\n");
        exit(1);
    }
    if (sc->mca == NULL){
        fprintf(stderr,"Must initialize markov chain approximation method\n");
        fprintf(stderr,"before initializing the dynamic program\n");
        exit(1);
    }
    if (sc->opt == NULL){
        fprintf(stderr,"Must add optimizer to c3sc before initializing\n");
        fprintf(stderr,"the dynamic program\n");
        exit(1);
    }
    
    assert (sc->type == IH);
    sc->dpih = dpih_alloc(beta,s,b);
    dpih_attach_mca(sc->dpih, sc->mca);
    dpih_attach_cost(sc->dpih, sc->cost);
    dpih_attach_policy(sc->dpih, sc->pol);
    dpih_attach_opt(sc->dpih,sc->opt);
}

void * c3sc_get_dp(struct C3SC * sc)
{
    assert (sc != NULL);
    assert (sc->type == IH);
    return sc->dpih;
}



/** \struct DPih
 *  \brief Infinite Horizon dynamic program
 *  \var DPih::mm
 *  markov chain approximation
 *  \var DPih::cost
 *  cost function
 *  \var DPih::pol
 *  policy
 *  \var DPih::beta
 *  >0 discount factor
 *  \var DPih::stagecost
 *  stage cost f(time,state,control,out,grad)
 *  \var DPih::boundcost
 *  boundary cost f(time,state,out)
 *  \var DPih::opt
 *  optimization options
 */
struct DPih
{
    struct MCA * mm;
    struct Cost * cost;
    struct Policy * pol;

    double beta; // discount factor
    int (*stagecost)(double,double*,double*,double*,double*);
    int (*boundcost)(double,double*,double*);

    struct c3Opt * opt;
};

///////////////////////////////////////////////////////////////

/**********************************************************//**
    Allocate an infinite horizon Dynamic program

    \param[in] beta - discount factor
    \param[in] s    - stagecost
    \param[in] b    - boundary cost

    \return dynamic program
**************************************************************/
struct DPih * 
dpih_alloc(double beta,
           int (*s)(double,double*,double*,double*,double*),
           int (*b)(double,double*,double*))
{
    struct DPih * dp = malloc(sizeof(struct DPih));
    if (dp == NULL){
        fprintf(stderr, "Allocating DPih failed\n");
        exit(1);
    }

    dp->mm = NULL;
    dp->cost = NULL;
    dp->pol = NULL;
    dp->beta = beta;
    dp->stagecost = s;
    dp->boundcost = b;

    dp->opt = NULL;
    return dp;
}

/**********************************************************//**
    Free dynamic program
**************************************************************/
void dpih_free(struct DPih * dp)
{
    if (dp != NULL){
        free(dp); dp = NULL;
    }
}

/**********************************************************//**
    Free dynamic program and associated markov model, cost,
    policy, and optimization algorithm
**************************************************************/
void dpih_free_deep(struct DPih * dp)
{
    if (dp != NULL){
        policy_free(dp->pol); dp->pol = NULL;
        mca_free_deep(dp->mm); dp->mm = NULL;
        cost_free(dp->cost); dp->cost = NULL;
        c3opt_free(dp->opt); dp->opt = NULL;
        dpih_free(dp); dp = NULL;
    }
}

/**********************************************************//**
    Attach a reference to a markov chain approximation
**************************************************************/
void dpih_attach_mca(struct DPih * dp, struct MCA * mm)
{
    assert(dp!= NULL);
    dp->mm = mm;
}

/**********************************************************//**
   Attach a reference to a cost function 
**************************************************************/
void dpih_attach_cost(struct DPih * dp, struct Cost * cost)
{
    assert(dp != NULL);
    dp->cost = cost;
}

/**********************************************************//**
   Attach a reference to a policy
**************************************************************/
void dpih_attach_policy(struct DPih * dp, struct Policy * pol)
{
    assert (dp != NULL);
    dp->pol = pol;
}

/**********************************************************//**
   Attach an optimization routine
**************************************************************/
void dpih_attach_opt(struct DPih * dp, struct c3Opt * opt)
{
    assert (dp != NULL);
    dp->opt = opt;
}

/**********************************************************//**
   Get the cost function
**************************************************************/
struct Cost * dpih_get_cost(struct DPih * dp)
{
    assert (dp != NULL);
    return dp->cost;
}

/**********************************************************//**
   Get the dynamics
**************************************************************/
struct Dyn * dpih_get_dyn(struct DPih * dp)
{
    assert (dp != NULL);
    return mca_get_dyn(dp->mm);
}

/**********************************************************//**
   Evaluate right hand side of Bellman equation at a given node *x*
   and for a given control *u*
**************************************************************/
double dpih_rhs(struct DPih * dp,double * x,double * u, double * grad)
{
    double dt;
    double * gdt = NULL;
    double t = 0.0;
    int absorb = 0;
    double val;
    size_t du;
    if (grad == NULL){
        val = mca_expectation(dp->mm,t,x,u,&dt,NULL,
                              cost_eval_bb,dp->cost,
                              &absorb,grad);
    }
    else{
        du = mca_get_du(dp->mm);
        gdt = calloc_double(du);
        val = mca_expectation(dp->mm,t,x,u,&dt,gdt,
                              cost_eval_bb,dp->cost,
                              &absorb,grad);
    }

    double out;
    if (absorb == 1){
        int res = dp->boundcost(t,x,&out);
        assert (res == 0);
    }
    else{
        double sc;

        if (grad == NULL){
            int res = dp->stagecost(t,x,u,&sc,NULL);
            assert (res == 0);
            out = exp(-dp->beta*dt)*val + dt*sc;
        }
        else{
            du = mca_get_du(dp->mm);
            double * gtemp = calloc_double(du);
            int res = dp->stagecost(t,x,u,&sc,gtemp);
            assert (res == 0);
            double ebt = exp(-dp->beta*dt);
            
            out = ebt*val + dt*sc;
            
            for (size_t jj = 0; jj < du; jj++){
                grad[jj] = -dp->beta*ebt*gdt[jj]*val +
                    ebt*grad[jj] + dt*gtemp[jj] + gdt[jj]*sc;
            }
            free(gtemp); gtemp = NULL;
        }
    }
    free(gdt); gdt = NULL;

    return out;
}

/**********************************************************//**
   Use the DP policy to compute the rhs of Bellman equation
**************************************************************/
double dpih_rhs_pol(struct DPih * dp, double * x)
{
    struct Control * u = NULL;
    double t = 0.0;
    int res = policy_eval(dp->pol,t,x,&u);
    assert (res == 0);

    double * uu = control_getu_ref(u);
    double val = dpih_rhs(dp,x,uu,NULL);
    
    control_free(u); u = NULL;
    return val;    
}

/**********************************************************//**
   Helper function for sub iteration of policy iteration
**************************************************************/
double dpih_rhs_bb(double * x, void * dp)
{
    double ret = dpih_rhs_pol(dp,x);
    return ret;
}

/**********************************************************//**
   Generate a new cost function by iterating Bellman equation
   with a *fixed* policy
**************************************************************/
struct Cost * dpih_iter_pol(struct DPih * dp,int verbose)
{
    struct Cost * oc = dp->cost;
    struct Cost * cost = cost_alloc(oc->d,oc->bds->lb,oc->bds->ub);
    cost_init_discrete(cost,oc->N,oc->x);
    cost_approx(cost,dpih_rhs_bb,dp,verbose);
    return cost;
}

struct DPX
{
    struct DPih * dp;
    double * x;
};
    
/**********************************************************//**
   Helper function for obtaining new cost by minimizing 
   Bellman Equation
**************************************************************/
double dpih_rhs_opt_bb(size_t du, double * u, double * grad, void * arg)
{
    (void)(du);
    struct DPX * dpx = arg;
    double val = dpih_rhs(dpx->dp,dpx->x,u,grad);
    return val;    
}

/**********************************************************//**
      Run optimizer and return optimal cost
**************************************************************/
double dpih_rhs_opt_cost(double * x,void * dp)
{

    struct DPX dpx;
    dpx.dp = dp;
    dpx.x = x;
    
    assert (dpx.dp->opt != NULL);

    size_t du = mca_get_du(dpx.dp->mm);
    double * ustart = calloc_double(du);
    
    double val = 0.0;
    struct c3Opt * opt = dpx.dp->opt;
    c3opt_add_objective(opt,dpih_rhs_opt_bb,&dpx);
    
    int res = c3opt_minimize(opt,ustart,&val);
    assert (res > -1);
     
    free(ustart); ustart = NULL;
    return val;
}

/**********************************************************//**
   Generate a new cost function by iterating Bellman equation
   with an *optimal* policy
**************************************************************/
struct Cost * dpih_iter_vi(struct DPih * dp,int verbose)
{
    struct Cost * oc = dp->cost;
    struct Cost * cost = cost_alloc(oc->d,oc->bds->lb,oc->bds->ub);
    cost_init_discrete(cost,oc->N,oc->x);
    cost_approx(cost,dpih_rhs_opt_cost,dp,verbose);
    return cost;
}

/**********************************************************//**
   For computing a policy
**************************************************************/
double dpih_rhs_opt_pol(double * x,size_t ind,void * dpin)
{
    struct DPX dpx;
    dpx.dp = dpin;
    dpx.x = x;
    
    assert (dpx.dp->opt != NULL);

    size_t du = mca_get_du(dpx.dp->mm);
    double * ustart = calloc_double(du);

    double val = 0.0;
    struct c3Opt * opt = dpx.dp->opt;
    c3opt_add_objective(opt,dpih_rhs_opt_bb,&dpx);
    
    int res = c3opt_minimize(opt,ustart,&val);
    assert (res > -1);
     
    val = ustart[ind];
    free(ustart); ustart = NULL;
    return val;
}


/**********************************************************//**
   Generate a new policy by iterating optimal Bellman equation
**************************************************************/
struct Policy * dpih_iter_vi_pol(struct DPih * dp,int verbose)
{
    assert (dp != NULL);
    assert (dp->cost != NULL);

    struct Cost * oc = dp->cost;
    double * lb = cost_get_lb(oc);
    double * ub = cost_get_ub(oc);
    size_t dx = mca_get_dx(dp->mm);
    size_t du = mca_get_du(dp->mm);

    struct Policy * pol = policy_alloc(dx,du);
    policy_set_bounds(pol,lb,ub);
    policy_init_discrete(pol,oc->N,oc->x);
    policy_approx(pol,dpih_rhs_opt_pol,dp,verbose);
    return pol;
}
