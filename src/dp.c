#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "c3.h"
#include "cdyn/src/simulate.h"

#include "util.h"
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
//    struct ImplictPolicy * pol;
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
//    sc->pol = NULL;
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
    Add bounds on state and initialize boundary conditions
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

    sc->bound = boundary_alloc(sc->dx,sc->lbx,sc->ubx);
}

/**********************************************************//**
    Set the external boundary type of a particular dimension
**************************************************************/
void c3sc_set_external_boundary(struct C3SC * sc, size_t dim,
                                char * type)
{
    boundary_external_set_type(sc->bound,dim,type);
}

/**********************************************************//**
    Add a rectangular obstacle
**************************************************************/
void c3sc_add_obstacle(struct C3SC * sc, double * center, double * lengths)
{
    boundary_add_obstacle(sc->bound,center,lengths);
    
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
                       int (*b)(double,const double*,const double*,
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
        /* printf("x[%zu] = ",ii);dprint(N[ii],x[ii]); */
        /* printf("(%G,%G)\n",sc->lbx[ii],sc->ubx[ii]); */
    }
    cost_init_discrete(sc->cost,N,x);

    size_t nobs = boundary_get_nobs(sc->bound);
    for (size_t jj = 0; jj < nobs; jj++){
        double * lb = boundary_obstacle_get_lb(sc->bound,jj);
        double * ub = boundary_obstacle_get_ub(sc->bound,jj);
        cost_add_nodes(sc->cost,lb,ub,3);
    }

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
    \param[in]     o    - obstacle cost (NULL if none)
**************************************************************/
void c3sc_init_dp(struct C3SC * sc, double beta,
                  int (*s)(double,double*,double*,double*,double*),
                  int (*b)(double,double*,double*),
                  int (*o)(double*,double*))
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
    sc->dpih = dpih_alloc(beta,s,b,o);
    dpih_attach_mca(sc->dpih, &(sc->mca));
    dpih_attach_cost(sc->dpih, &(sc->cost));
    dpih_attach_opt(sc->dpih, &(sc->opt));
}


/**********************************************************//**
    Load a cost function
**************************************************************/
int c3sc_cost_load(struct C3SC * sc, char * filename )
{
    assert (sc != NULL);
    assert (sc->cost != NULL);
    int load_success = cost_load(sc->cost,filename);
    return load_success;
}

/**********************************************************//**
    Approximate a cost function
**************************************************************/
int c3sc_cost_approx(struct C3SC * sc, double (*f)(double *, void *),
                     void * args, int verbose,
                     const struct ApproxArgs * aargs)
{
    assert (sc != NULL);
    assert (sc->cost != NULL);
    cost_approx(sc->cost,f,args,verbose,aargs);
    return 0;
}


/**********************************************************//**
    Update the cost function to be that corresponding 
    to the currently stored implicit policy
**************************************************************/
void c3sc_pol_solve(struct C3SC * sc, size_t max_solve_iter,
                    double solve_tol,
                    int verbose,
                    const struct ApproxArgs * aargs)
{
    struct ImplicitPolicy * pol = c3sc_create_implicit_policy(sc);
    dpih_iter_pol_solve(sc->dpih,pol,max_solve_iter,solve_tol,verbose,aargs);
    sc->cost = dpih_get_cost(sc->dpih);
    implicit_policy_free(pol);
}

/**********************************************************//**
   Take a step of value iteration        
**************************************************************/
double c3sc_iter_vi(struct C3SC * sc,int verbose,const struct ApproxArgs * aargs,
                    struct C3SCDiagnostic * diag)
{
    struct Cost * newcost = dpih_iter_vi(sc->dpih,verbose-1,aargs,diag);
    struct Cost * oldcost = dpih_get_cost(sc->dpih);
    double diff = cost_norm2_diff(newcost,oldcost);
    dpih_attach_cost_ow(sc->dpih,newcost);
    sc->cost = dpih_get_cost(sc->dpih);
    cost_free(newcost); newcost = NULL;
    return diff;
}

/**********************************************************//**
    Get a reference to the dynamic programming problem
**************************************************************/
void * c3sc_get_dp(const struct C3SC * sc)
{
    assert (sc != NULL);
    assert (sc->type == IH);
    return sc->dpih;
}

/**********************************************************//**
    Get a reference to the cost
**************************************************************/
struct Cost * c3sc_get_cost(const struct C3SC * sc)
{
    assert (sc != NULL);
    return sc->cost;
}

/**********************************************************//**
    Get the control size
**************************************************************/
size_t c3sc_get_du(const struct C3SC * sc)
{
    assert (sc != NULL);
    return sc->du;
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
 *  \var DPih::obscost
 *  cost of hitting obstacle f(state,out)
 *  \var DPih::opt
 *  optimization options
 */
struct DPih
{
    struct MCA * mm;
    struct Cost * cost;
    //  struct Policy * pol;

    double beta; // discount factor
    int (*stagecost)(double,double*,double*,double*,double*);
    int (*boundcost)(double,double*,double*);
    int (*obscost)(double*,double*);

    struct c3Opt * opt;
};

///////////////////////////////////////////////////////////////

/**********************************************************//**
    Allocate an infinite horizon Dynamic program

    \param[in] beta - discount factor
    \param[in] s    - stagecost
    \param[in] b    - boundary cost
    \param[in] obs  - cost of being in obstacle

    \return dynamic program
**************************************************************/
struct DPih * 
dpih_alloc(double beta,
           int (*s)(double,double*,double*,double*,double*),
           int (*b)(double,double*,double*),
           int (*obs)(double *,double*))
{
    struct DPih * dp = malloc(sizeof(struct DPih));
    if (dp == NULL){
        fprintf(stderr, "Allocating DPih failed\n");
        exit(1);
    }

    dp->mm = NULL;
    dp->cost = NULL;
//    dp->pol = NULL;
    dp->beta = beta;
    dp->stagecost = s;
    dp->boundcost = b;
    dp->obscost = obs;

    dp->opt = NULL;
    return dp;
}

/**********************************************************//**
    Interpolate a dynamic program
**************************************************************/
struct DPih * dpih_interp2(struct DPih * dp, int inhalf)
{
    struct DPih * newdp = dpih_alloc(dp->beta,dp->stagecost,
                                     dp->boundcost,dp->obscost);
    newdp->opt = c3opt_copy(dp->opt);

    newdp->cost = cost_copy_deep(dp->cost);
    cost_interp_inhalf(newdp->cost,inhalf);
    size_t d = cost_get_d(dp->cost); 
   
    double * newh = calloc_double(d);
    cost_get_h(newdp->cost,newh);
    
    newdp->mm = mca_copy_deep(dp->mm);
    mca_set_newh(newdp->mm,newh);
        
    free(newh); newh = NULL;

    return newdp;
}

/**********************************************************//**
    Copy dynamic program
**************************************************************/
struct DPih * dpih_copy_deep(struct DPih * dp)
{
    if (dp == NULL){
        return NULL;
    }

    struct DPih * newdp = dpih_alloc(dp->beta,dp->stagecost,
                                     dp->boundcost,dp->obscost);

    newdp->mm = mca_copy_deep(dp->mm);
    newdp->cost = cost_copy_deep(dp->cost);
    newdp->opt = c3opt_copy(dp->opt);
    
    return newdp;
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
//        policy_free(dp->pol); dp->pol = NULL;
        mca_free_deep(dp->mm); dp->mm = NULL;
        cost_free(dp->cost); dp->cost = NULL;
        c3opt_free(dp->opt); dp->opt = NULL;
        dpih_free(dp); dp = NULL;
    }
}

/**********************************************************//**
    Attach a reference to a markov chain approximation
**************************************************************/
void dpih_attach_mca(struct DPih * dp, struct MCA ** mm)
{
    assert(dp!= NULL);
    dp->mm = *mm;
}

/**********************************************************//**
   Attach a reference to a cost function 
**************************************************************/
void dpih_attach_cost(struct DPih * dp, struct Cost ** cost)
{
    assert(dp != NULL);
    dp->cost = *cost;
}

/**********************************************************//**
   Overwrite a previous cost function
**************************************************************/
void dpih_attach_cost_ow(struct DPih * dp, struct Cost * cost)
{
    assert(dp != NULL);
    cost_free(dp->cost); dp->cost = NULL;
    dp->cost = cost_copy_deep(cost);
}

/**********************************************************//**
   Attach an optimization routine
**************************************************************/
void dpih_attach_opt(struct DPih * dp, struct c3Opt ** opt)
{
    assert (dp != NULL);
    dp->opt = *opt;
}

/**********************************************************//**
   Get the state dimension
**************************************************************/
size_t dpih_get_d(const struct DPih * dp)
{
    assert (dp != NULL);
    return mca_get_dx(dp->mm);
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
    size_t du = mca_get_du(dp->mm);
    size_t dx = mca_get_dx(dp->mm);
    /* printf("eval rhs\n"); */
    /* printf("x = "); dprint(mca_get_dx(dp->mm),x); */
    /* printf("u = "); dprint(du,u); */
    double dt;
    double * gdt = NULL;
    double t = 0.0;
    struct BoundInfo * bi = NULL;
    double val;

    /* printf("get expectation \n"); */

    int info;
    if (grad == NULL){
        val = mca_expectation(dp->mm,t,x,u,&dt,NULL,
                              cost_eval_bb,dp->cost,
                              &bi,grad,&info);
    }
    else{
        gdt = calloc_double(du);
        val = mca_expectation(dp->mm,t,x,u,&dt,gdt,
                              cost_eval_bb,dp->cost,
                              &bi,grad,&info);
    }
    if (info != 0){
        fprintf(stderr, "Warning: mca_expectation call from dpih_rhs");
        fprintf(stderr, " resulted in an error\n");
        fprintf(stderr, "\tx = ");
        for (size_t ii = 0; ii < dx; ii++ ){
            fprintf(stderr,"%G ",x[ii]);
        }
        fprintf(stderr,"\n\tControl = ");
        for (size_t ii = 0; ii < du; ii++){
            fprintf(stderr,"%G ",u[ii]);
        }
        fprintf(stderr,"\n");
        fprintf(stderr,"\tBoundary info is ?????\n");
        exit(1);
    }
    
    /* printf("got expectation\n"); */
    double out;
    int absorb = bound_info_absorb(bi);
    if (absorb == 1){
//        printf("MASSIVE ERROR!\n");
        int in_obstacle = bound_info_get_in_obstacle(bi);
        if (in_obstacle < 0){
            int res = dp->boundcost(t,x,&out);
            assert (res == 0);
        }
        else{
            assert (dp->obscost != NULL);
            int res = dp->obscost(x,&out);
            assert(res == 0);
        }
        if (grad != NULL){
            for (size_t ii = 0; ii < du; ii++ ){
                grad[ii] = 0.0;
            }
        }
        /* printf("\n\nx = ");dprint(2,x); */
        /* printf("in obstacle %d\n",in_obstacle); */
        /* printf("out = %G\n",out); */

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
    bound_info_free(bi); bi = NULL;
    /* printf("evaluated it\n"); */
    return out;
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
    /* if (x[2] < -1.6){ */
    /*     if (x[2] > -2.0){ */
    /*         printf("shouldn't be calling this\n"); */
    /*         dprint(3,x); */
    /*         exit(1); */
    /*     } */
    /* } */
    struct DPX dpx;
    dpx.dp = dp;
    dpx.x = x;
    
    assert (dpx.dp->opt != NULL);

    size_t du = mca_get_du(dpx.dp->mm);
    double * ustart = calloc_double(du);
    /* double out = dpih_rhs(dpx.dp,x,ustart,NULL); */
    /* free(ustart); */
    /* return out; */
    
    double val = 0.0;
    struct c3Opt * opt = dpx.dp->opt;
    c3opt_add_objective(opt,dpih_rhs_opt_bb,&dpx);
    size_t dx = mca_get_dx(dpx.dp->mm);
//    printf("before = ");
//    dprint(dx,x);
    int res = c3opt_minimize(opt,ustart,&val);
    if (res < -1){
        printf("max iter reached in optimization res=%d\n",res);

        printf("x = ");
        dprint(dx,x);
        dprint(du,ustart);
        for (size_t ii = 0; ii < du; ii++){
            ustart[ii] = 0.001;
        }
        val = 0.0;
        printf("restart with verbose\n");
        c3opt_set_verbose(opt,2);
        int res2 = c3opt_minimize(opt,ustart,&val);
        printf("res2 = %d\n",res2);
    }
    // assert (res > -1);

    /* if (ustart[0] > 1e-2){ */
    /*     printf("x = "); dprint(2,x); */
    /*     printf("u = %3.15G\n",ustart[0]); */
    /* } */

    /* printf("x = "); dprint(dx, x); */
    /* printf("optimal u = %G\n",ustart[0]); */
    free(ustart); ustart = NULL;
    return val;
}

struct DPPOL
{
    struct DPih * dp;
    struct ImplicitPolicy * pol;
    double weight;
};
    
/**********************************************************//**
      Bellman RHS with fixed policy
**************************************************************/
double dpih_rhs_pol_cost(double * x,void * dp)
{

    struct DPPOL * dppol = dp;
    struct DPX dpx;
    dpx.dp = dppol->dp;
    dpx.x = x;

    size_t du = mca_get_du(dpx.dp->mm);
    double * u = calloc_double(du);
    int res = implicit_policy_eval(dppol->pol,0.0,x,u);
    assert (res == 0);
    double val = dpih_rhs(dppol->dp,x,u,NULL);
    free(u); u = NULL;
    return val;
}

/**********************************************************//**
      Bellman RHS with fixed policy
**************************************************************/
double dpih_rhs_pol_cost_weight(double * x,void * dp)
{

    struct DPPOL * dppol = dp;
    struct DPX dpx;
    dpx.dp = dppol->dp;
    dpx.x = x;

    size_t du = mca_get_du(dpx.dp->mm);
    double * u = calloc_double(du);
    int res = implicit_policy_eval(dppol->pol,0.0,x,u);
    assert (res == 0);
    double val = dpih_rhs(dppol->dp,x,u,NULL);
    double oldval = 0;
    res = cost_eval(dpx.dp->cost,0,x,&oldval);
    assert (res == 0);
    val = dppol->weight * val + (1.0 - dppol->weight)*oldval;
    free(u); u = NULL;

    return val;
}

/**********************************************************//**
   Generate a new cost function by iterating Bellman equation
   with an *optimal* policy
**************************************************************/
struct Cost * dpih_iter_vi(struct DPih * dp,int verbose,
                           const struct ApproxArgs * aargs,
                           struct C3SCDiagnostic * diag)
{

//    printf("do iter vi\n");
    struct Cost * cost = cost_copy_deep(dp->cost);
//    printf("done\n");
//    struct Cost * cost = cost_alloc(d,lb,ub);
//    cost_init_discrete(cost,oc->N,oc->x);

    size_t d = cost_get_d(cost);
    struct FunctionMonitor * fm = NULL;
    fm = function_monitor_initnd(dpih_rhs_opt_cost,dp,d,1000*d);
    cost_approx(cost,function_monitor_eval,fm,verbose,aargs);
    
    if (diag != NULL)
    {
        c3sc_diagnostic_vi_update(diag,cost,dp,fm);
    }
    
    function_monitor_free(fm); fm = NULL;

//    cost_approx(cost,dpih_rhs_opt_cost,dp,verbose,cross_tol,round_tol,kickrank);
    return cost;
}

/**********************************************************//**
   Generate a new cost function by iterating Bellman equation
   with a fixed policy
**************************************************************/
struct Cost * dpih_iter_pol(struct DPih * dp, struct ImplicitPolicy * pol,
                            int verbose, const struct ApproxArgs * aargs)
{
    struct DPPOL dppol;
    dppol.dp = dp;
    dppol.pol = pol;
    struct Cost * oc = dp->cost;

    /* size_t d = cost_get_d(oc); */
    /* double * lb = cost_get_lb(oc); */
    /* double * ub = cost_get_ub(oc); */
    /* struct Cost * cost = cost_alloc(d,lb,ub); */
    /* cost_init_discrete(cost,oc->N,oc->x); */

    struct Cost * cost = cost_copy_deep(oc);
//    printf("got first one!\n");
    cost_approx(cost,dpih_rhs_pol_cost,&dppol,verbose,aargs);
    return cost;
}

/**********************************************************//**
   Generate a new cost function by iterating a weighted
   version of the Bellman equation with a fixed policy
**************************************************************/
struct Cost * 
dpih_iter_pol_weight(struct DPih * dp, struct ImplicitPolicy * pol, 
                     double weight,
                     int verbose, const struct ApproxArgs * aargs)
{
    struct DPPOL dppol;
    dppol.dp = dp;
    dppol.pol = pol;
    dppol.weight = weight;
    struct Cost * oc = dp->cost;

    /* size_t d = cost_get_d(oc); */
    /* double * lb = cost_get_lb(oc); */
    /* double * ub = cost_get_ub(oc); */
    /* struct Cost * cost = cost_alloc(d,lb,ub); */
    /* cost_init_discrete(cost,oc->N,oc->x); */

    struct Cost * cost = cost_copy_deep(oc);
//    printf("got first one!\n");
    cost_approx(cost,dpih_rhs_pol_cost_weight,&dppol,verbose,aargs);
    return cost;
}

/**********************************************************//**
   Solve for the cost of a particular policy  
**************************************************************/
void
dpih_iter_pol_solve(struct DPih * dp,
                    struct ImplicitPolicy * pol,
                    size_t max_solve_iter, double solve_tol,
                    int verbose,
                    const struct ApproxArgs * aargs)
{
    double normprev = 0.0;

    struct Cost * tcost = NULL;
    double prevrat = 0;
    for (size_t jj = 0; jj < max_solve_iter; jj++){

        tcost = dpih_iter_pol(dp,pol,verbose-1,aargs);

        double normval = cost_norm2(tcost);
        double rat = fabs(normprev/normval);
        struct Cost * cost = dpih_get_cost(dp);
        double diff = cost_norm2_diff(tcost,cost);
        if (verbose > 0){
            if ((jj+1) % 1 == 0){
                printf("\t jj=%zu, normval=%G, ratio=%G\n",jj,diff/normval,rat);
            }
        }
        if (jj > 1){
            if ((prevrat < 1) && (rat > 1)){
                break;
            }
            if ((prevrat > 1) && (rat < 1)){
                break;
            }
        }
        prevrat = rat;
        normprev = normval;
        dpih_attach_cost_ow(dp,tcost);
        cost_free(tcost); tcost = NULL;
        if (diff/normval <  solve_tol){
            break;
        }
    }
    
//    return tcost;
}


/** \struct ImplicitPolicy
 *  \brief Stores information needed to evaluate implicit control
 *  \var ImplicitPolicy::dpalloc
 *  flag to indicate whether need to free dp on cleanup
 *  \var ImplicitPolicy::dp
 *  infinite horizon dp problem
 *  \var ImplicitPolicy::du
 *  size of control
 *  \var ImplicitPolicy::fm
 *  store previously calculated controls
 *  \var ImplicitPolicy::transform
 *  transformation function in case need to perform state transformation before evaluating policy
 */
struct ImplicitPolicy
{
    int dpalloc;
    struct DPih * dp;
    size_t du;
    struct HashtableCpair ** fm;

    size_t dx;
    void (*transform)(size_t,const double*,double*);
};

/**********************************************************//**
   Allocate memory for policy
**************************************************************/
struct ImplicitPolicy * implicit_policy_alloc()
{
    struct ImplicitPolicy * ip = malloc(sizeof(struct ImplicitPolicy));
    if (ip == NULL){
        fprintf(stderr, "Failure allocating memory for implicit policy\n");
        exit(1);
    }
    ip->dpalloc = 0;
    ip->dp = NULL;
    ip->fm = NULL;

    ip->dx = 0;
    ip->transform = NULL;
    return ip;
}

/**********************************************************//**
    Create an implicit policy
**************************************************************/
struct ImplicitPolicy * c3sc_create_implicit_policy(struct C3SC * sc)
{
    assert (sc != NULL);
    assert (sc->type == IH);
    struct ImplicitPolicy * ip = implicit_policy_alloc();
    struct DPih * dp = dpih_copy_deep(sc->dpih);
    implicit_policy_set_dp(ip,dp);


    ip->du = c3sc_get_du(sc);
    ip->fm = malloc(ip->du * sizeof(struct HashtableCpair *));
    if (ip->fm == NULL){
        fprintf(stderr,"Cannot allocate memory for implicit policy\n");
        exit(1);
    }
    size_t d = dpih_get_d(dp);
    for (size_t ii = 0; ii < ip->du; ii++)
    {
        ip->fm[ii] = create_hashtable_cp(d*1000);
    }
    /* printf("created implicit policy!\n"); */

    return ip;
}

/**********************************************************//**
    Add transform to the policy
**************************************************************/
void implicit_policy_add_transform(struct ImplicitPolicy * pol,size_t dx, void (*f)(size_t, const double *, double *))
{
    pol->dx = dx;
    pol->transform = f;
}


/**********************************************************//**
   Free policy
**************************************************************/
void implicit_policy_free(struct ImplicitPolicy * ip)
{
    if (ip != NULL){
        if (ip->dpalloc == 1){
            dpih_free_deep(ip->dp); ip->dp = NULL;
        }
        for (size_t ii = 0; ii < ip->du; ii++){
            free_hashtable_cp(ip->fm[ii]); ip->fm[ii] = 0;
        }
        free(ip->fm); ip->fm = NULL;
        free(ip); ip = NULL;
    }
}

/**********************************************************//**
   SetDP
**************************************************************/
void implicit_policy_set_dp(struct ImplicitPolicy * ip, struct DPih * dp)
{
    if (ip->dpalloc == 1){
        if (ip->dp != NULL){
            dpih_free_deep(ip->dp); ip->dp = NULL;
        }
    }
    ip->dpalloc = 1;
    ip->dp = dp;
}

/**********************************************************//**
   Evaluate a policy
   
   \param[in]     ip  - policy to evaluate
   \param[in]     t   - time at which to evaluate
   \param[in]     xin - location at which to evaluate
   \param[in,out] u   - resulting control

   \return 0 - success, else failure
**************************************************************/
int implicit_policy_eval(struct ImplicitPolicy * ip,double t,
                         const double * xin, double * u)
{
    (void)(t);

    /* printf("here\n"); */
    struct DPX dpx;
    dpx.dp = ip->dp;
    if (ip->transform == NULL){
        dpx.x = (double *) xin;
    }
    else{
        dpx.x = calloc_double(ip->dx);
        ip->transform(ip->dx,xin,dpx.x);
        /* printf("transformed!\n"); */
        /* printf("ip->dx = "); */
        /* dprint(ip->dx,dpx.x); */
    }


    size_t dx = mca_get_dx(dpx.dp->mm);
    size_t du = mca_get_du(dpx.dp->mm);

    char * ser = serialize_darray_to_text(dx,dpx.x);
    char * sval = lookup_key(ip->fm[0],ser);
    /* if (1 == 0){ */
    if (sval != NULL){
        /* printf("sval exists\n"); */
        u[0] = deserialize_double_from_text(sval);
        for (size_t ii = 1; ii < du; ii++){
            char * sval2 = lookup_key(ip->fm[ii],ser);
            u[ii] = deserialize_double_from_text(sval2);
            free(sval2); sval2 = NULL;
        }
    }
    else{
        /* printf("sval doesnt exists\n"); */
        assert (dpx.dp->opt != NULL);
        assert (dpx.dp->cost != NULL);
        assert (dpx.dp->mm != NULL);
        double val = 0.0;
        struct c3Opt * opt = dpx.dp->opt;
        c3opt_add_objective(opt,dpih_rhs_opt_bb,&dpx);
        /* printf("minimize \n"); */
        int res = c3opt_minimize(opt,u,&val);
        /* printf("minimized!\n"); */
        if (res < -1){
            printf("max iter reached in optimization res=%d\n",res);

            printf("x = ");
            dprint(dx,dpx.x);
            dprint(du,u);
            for (size_t ii = 0; ii < du; ii++){
                u[ii] = 0.0;
            }
            val = 0.0;
            printf("restart with verbose\n");
            c3opt_set_verbose(opt,2);
            int res2 = c3opt_minimize(opt,u,&val);
            printf("res2 = %d\n",res2);
        }

        for (size_t ii = 0; ii < du; ii++){
            char * sval2 = serialize_double_to_text(u[ii]);
            struct Cpair * cp = cpair_create(ser,sval2);
            add_cpair(ip->fm[ii],cp);
            cpair_free(cp); cp = NULL;
            free(sval2);
        }
    }
    free(ser); ser = NULL;
    free(sval); sval = NULL;

//    printf("we are done!\n");
    if (ip->transform != NULL){
        free(dpx.x); dpx.x = NULL;
    }
    int res_pol = 0;
    return res_pol;
}


/**********************************************************//**
   Interface for integrator
**************************************************************/
int implicit_policy_controller(double t,const double * x, double * u, void * args)
{

    struct ImplicitPolicy * ip = args;
    return implicit_policy_eval(ip,t,x,u);
}



/////////////////////////////////////////////////////////////////////////
// diagnostic stuff
////////////////////////////////////////////////////////////////////////
struct C3SCDiagnostic
{
    size_t itervi;
    size_t iterpi;
    struct Trajectory * traj_vi;
};

struct C3SCDiagnostic * c3sc_diagnostic_init()
{
    struct C3SCDiagnostic * diag = malloc(sizeof(struct C3SCDiagnostic));
    if (diag == NULL){
        fprintf(stderr,"Failure allocating diagnostic arguments\n");
        exit(1);
    }
    diag->itervi = 0;
    diag->iterpi = 0;
    diag->traj_vi = NULL;

    return diag;
}

void c3sc_diagnostic_free(struct C3SCDiagnostic * diag)
{
    if (diag != NULL){
        trajectory_free(diag->traj_vi); diag->traj_vi = NULL;
        free(diag); diag = NULL;
    }
}

int c3sc_diagnostic_save(struct C3SCDiagnostic * diag, 
                         char * filename, size_t nprec)
{
    FILE * fp = fopen(filename,"w");
    if (fp == NULL){
        return 1;
    }
    trajectory_print(diag->traj_vi,fp,nprec);
    fclose(fp);
    return 0;
}

void c3sc_diagnostic_vi_update(struct C3SCDiagnostic * diag,
                               struct Cost * newcost,
                               struct DPih * dp,
                               struct FunctionMonitor * fm)
{
    struct Cost * old = dpih_get_cost(dp);
    size_t d = cost_get_d(old);

    double norm = cost_norm2(newcost);
    double diff = cost_norm2_diff(old,newcost);
    size_t * ranks = cost_get_ranks(newcost);
    size_t ntot = cost_get_size(newcost);
    size_t nevals = nstored_hashtable_cp(fm->evals);
    double frac = (double) nevals / ntot;

    size_t ndim = d+1+3;
    double * pt = calloc_double(ndim);
    pt[0] = norm;
    pt[1] = diff;
    pt[2] = frac;
    for (size_t ii = 0; ii < d+1; ii++){
        pt[ii+3] = (double) ranks[ii];
    }
    
    diag->itervi += 1;
    trajectory_add(&(diag->traj_vi),ndim,0,(double)(diag->itervi),pt,NULL);

    free(pt); pt = NULL;

}

