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
    int (*stagecost)(double,const double*,const double*,double*,double*);
    int (*boundcost)(double,const double*,double*);
    int (*obscost)(const double*,double*);

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
           int (*s)(double,const double*,const double*,double*,double*),
           int (*b)(double,const double*,double*),
           int (*obs)(const double *,double*))
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
/* struct DPih * dpih_interp2(struct DPih * dp, int inhalf) */
/* { */
/*     struct DPih * newdp = dpih_alloc(dp->beta,dp->stagecost, */
/*                                      dp->boundcost,dp->obscost); */
/*     newdp->opt = c3opt_copy(dp->opt); */

/*     newdp->cost = cost_copy_deep(dp->cost); */
/*     cost_interp_inhalf(newdp->cost,inhalf); */
/*     size_t d = cost_get_d(dp->cost);  */
   
/*     double * newh = calloc_double(d); */
/*     cost_get_h(newdp->cost,newh); */
    
/*     newdp->mm = mca_copy_deep(dp->mm); */
/*     mca_set_newh(newdp->mm,newh); */
        
/*     free(newh); newh = NULL; */

/*     return newdp; */
/* } */

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
   Get the MCA
**************************************************************/
struct MCA * dpih_get_mca(struct DPih * dp)
{
    assert(dp!= NULL);
    return dp->mm;
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
   Get the optimization struct
**************************************************************/
struct c3Opt * dpih_get_opt(struct DPih * dp)
{
    assert (dp != NULL);
    return dp->opt;
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


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
////////////////////    Solvers     //////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



/**********************************************************//**
   Evaluate right hand side of Bellman equation at a given node *x*
   and for a given control *u*
**************************************************************/
double dpih_bellman(struct DPih * dp,
                    struct Node * node,
                    struct BoundInfo * bi,
                    const double * u,
                    double * grad)
{
    
    size_t du = mca_get_du(dp->mm);
    size_t dx = mca_get_dx(dp->mm);
    double t = 0.0;

    double dt;
    double * gdt = NULL;

    double val;

    /* printf("rhs bellman \n"); */
    int info;
    if (grad == NULL){
        val = mca_expectation_cost2(dp->mm,t,node,u,&dt,NULL,
                                    grad,&info);
    }
    else{
        gdt = calloc_double(du);
        val = mca_expectation_cost2(dp->mm,t,node,u,&dt,gdt,
                                    grad,&info);
    }
    /* printf("got it\n"); */
    if (info != 0){
        fprintf(stderr, "Warning: mca_expectation call from dpih_rhs");
        fprintf(stderr, " resulted in an error\n");
        fprintf(stderr, "\tx = ");
        for (size_t ii = 0; ii < dx; ii++ ){
            fprintf(stderr,"%G ",node->x[ii]);
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
            int res = dp->boundcost(t,node->x,&out);
            assert (res == 0);
        }
        else{
            assert (dp->obscost != NULL);
            int res = dp->obscost(node->x,&out);
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
            int res = dp->stagecost(t,node->x,u,&sc,NULL);
            assert (res == 0);
            out = exp(-dp->beta*dt)*val + dt*sc;
        }
        else{
            du = mca_get_du(dp->mm);
            double * gtemp = calloc_double(du);
            int res = dp->stagecost(t,node->x,u,&sc,gtemp);
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
    /* bound_info_free(bi); bi = NULL; */
    /* printf("evaluated it\n"); */
    return out;
}

/**********************************************************//**
      Bellman RHS with fixed policy
**************************************************************/
double dpih_bellman_policy(struct DPih * dp,
                           struct Node * node,
                           struct BoundInfo * bi,
                           struct ImplicitPolicy * pol)
{

    size_t du = mca_get_du(dp->mm);
    double * u = calloc_double(du);
    int res = implicit_policy_eval(pol,0.0,node->x,u);
    assert (res == 0);
    double val = dpih_bellman(dp,node,bi,u,NULL);
    free(u); u = NULL;
    return val;
}

/**********************************************************//**
   Helper function for obtaining new cost by minimizing 
   Bellman Equation
**************************************************************/
double dpih_bellman_evalu(size_t du, double * u, double * grad, void * arg)
{
    (void)(du);
    struct DPX2 * dpx = arg;
    double val = dpih_bellman(dpx->dp,dpx->node,dpx->bi,u,grad);
    return val;    
}

/**********************************************************//**
      Run optimizer and return optimal cost
**************************************************************/
double dpih_bellman_minimize(const double * x,void * dp)
{

    struct DPX2 dpx;
    dpx.dp = dp;
    assert (dpx.dp->opt != NULL);

    /* printf("start\n"); */
    size_t dx = mca_get_dx(dpx.dp->mm);
    size_t du = mca_get_du(dpx.dp->mm);
    /* printf("x = "); dprint(dx,x); */
    struct Boundary * bound = mca_get_boundary(dpx.dp->mm);
    double * h = mca_get_h(dpx.dp->mm);
    double time = 0.0;
    
    struct BoundInfo * bi = boundary_type(bound,time,x);
    /* printf("\n\n\n\n init node dx=%zu\n",dx); dprint(dx,x); */
    struct Node * node = node_init(dx,x,h,dpx.dp->cost,bi);
    /* printf("node->d = %zu\n",node->d); */
    /* printf("node->use = ");  iprint(2*dx,node->use); */
    /* printf("node->c = %G\n",node->c); */
    /* printf("node->pert = "); dprint(2*dx,node->pert); */
    /* printf("node->cost = "); dprint(2*dx,node->cost); */
    dpx.node = node;
    dpx.bi = bi;

    struct c3Opt * opt = dpx.dp->opt;    
    c3opt_add_objective(opt,dpih_bellman_evalu,&dpx);

    double * ustart = calloc_double(du);
    double * utemp = calloc_double(du);
    double * lb = c3opt_get_lb(opt);
    double * ub = c3opt_get_ub(opt);

    size_t nsamples = 20;
    for (size_t ii = 0; ii < du; ii++){
        utemp[ii] = (ub[ii]-lb[ii])*randu() + lb[ii];
    }
    double mincost = c3opt_eval(opt,utemp,NULL);;
    memmove(ustart,utemp,du*sizeof(double));
    /* printf("(cost = %G for control \n \t",mincost); */
    /* dprint(du,utemp); */
    for (size_t kk = 0; kk < nsamples; kk++){
      
        for (size_t ii = 0; ii < du; ii++){
            utemp[ii] = (ub[ii]-lb[ii])*randu() + lb[ii];
        }
        double out = c3opt_eval(opt,utemp,NULL);
        /* printf("(cost = %G for control \n \t",out); */
        /* dprint(du,utemp); */
        if (out < mincost){
            mincost = out;
            memmove(ustart,utemp,du*sizeof(double));
        }

    }

    double val = mincost;
    int res = 1;
    /* double val = 0.0; */
    /* printf("ustart = "); */
    /* dprint(du,ustart); */
    /* int res = c3opt_minimize(opt,ustart,&val); */
    /* printf("minimized res = %d,val = %G,mincost=%G\n",res,val,mincost); */
    /* dprint(du,ustart); */
    /* assert (val <= mincost); */
    free(utemp); utemp = NULL;
    /* exit(1); */


//    printf("before = ");
//    dprint(dx,x);

    /* printf("do minimize\n"); */

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
    free(ustart); ustart = NULL;
    
    node_free(node); node = NULL;
    bound_info_free(bi); bi = NULL;
    /* printf("returnign minimize val = %G\n",val); */
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


    struct Cost * cost = cost_copy_deep(dp->cost);

    size_t d = cost_get_d(cost);
    struct FunctionMonitor * fm = NULL;
    fm = function_monitor_initnd(dpih_bellman_minimize,dp,d,1000*d);
    /* printf("approx cost\n"); */
    cost_approx(cost,function_monitor_eval,fm,verbose,aargs);
    /* printf("costed!\n"); */
    
    if (diag != NULL)
    {
        c3sc_diagnostic_vi_update(diag,cost,dp,fm);
    }
    
    function_monitor_free(fm); fm = NULL;

    return cost;
}

/**********************************************************//**
      Run optimizer and return optimal cost
**************************************************************/
double dpih_bellman_policy_iterate(const double * x,void * dp)
{

    struct DPPOL * dpx = dp;
    
    size_t dx = mca_get_dx(dpx->dp->mm);
    /* size_t du = mca_get_du(dpx->dp->mm); */
    struct Boundary * bound = mca_get_boundary(dpx->dp->mm);
    double * h = mca_get_h(dpx->dp->mm);
    double time = 0.0;
    
    struct BoundInfo * bi = boundary_type(bound,time,x);
    struct Node * node = node_init(dx,x,h,dpx->dp->cost,bi);

    double val = dpih_bellman_policy(dpx->dp,node,bi,dpx->pol);
    node_free(node); node = NULL;
    bound_info_free(bi); bi = NULL;
    return val;
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
    struct Cost * cost = cost_copy_deep(oc);
    
    cost_approx(cost,dpih_bellman_policy_iterate,&dppol,verbose,aargs);
    
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
    /* double prevrat = 0; */
    for (size_t jj = 0; jj < max_solve_iter; jj++){

        tcost = dpih_iter_pol(dp,pol,verbose-1,aargs);

        double normval = cost_norm2(tcost);
        double rat = fabs(normprev/normval);
        struct Cost * cost = dpih_get_cost(dp);
        double diff = cost_norm2_diff(tcost,cost);
        if (verbose > 0){
            if ((jj+1) % 1 == 0){
                printf("\t POLICY ITERATION (%zu\\%zu):\n",jj+1,max_solve_iter);
                printf("\t \t L2 Difference between iterates    = %3.5E\n ",diff);
                printf("\t \t L2 Norm of current value function = %3.5E\n", normval);
                printf("\t \t Relative L2 Cauchy difference     = %3.5E\n", diff/normval);
                printf("\t \t Ratio: L2 Norm Prev / L2 Norm Cur = %G\n",rat);
            }
        }
        if (jj > 1){
            /* if ((prevrat < 1) && (rat > 1)){ */
            /*     break; */
            /* } */
            /* if ((prevrat > 1) && (rat < 1)){ */
            /*     break; */
            /* } */
        }
        /* prevrat = rat; */
        normprev = normval;
        dpih_attach_cost_ow(dp,tcost);
        cost_free(tcost); tcost = NULL;
        if (diff/normval <  solve_tol){
            break;
        }
    }
    
//    return tcost;
}

/* /\**********************************************************\//\** */
/*       Bellman RHS with fixed policy */
/* **************************************************************\/ */
/* double dpih_rhs_pol_cost_weight(const double * x,void * dp) */
/* { */

/*     struct DPPOL * dppol = dp; */
/*     struct DPX dpx; */
/*     dpx.dp = dppol->dp; */
/*     dpx.x = x; */

/*     size_t du = mca_get_du(dpx.dp->mm); */
/*     double * u = calloc_double(du); */
/*     int res = implicit_policy_eval(dppol->pol,0.0,x,u); */
/*     assert (res == 0); */
/*     double val = dpih_rhs(dppol->dp,x,u,NULL); */
/*     double oldval = 0; */
/*     res = cost_eval(dpx.dp->cost,0,x,&oldval); */
/*     assert (res == 0); */
/*     val = dppol->weight * val + (1.0 - dppol->weight)*oldval; */
/*     free(u); u = NULL; */

/*     return val; */
/* } */

/**********************************************************//**
   Generate a new cost function by iterating a weighted
   version of the Bellman equation with a fixed policy
**************************************************************/
/* struct Cost *  */
/* dpih_iter_pol_weight(struct DPih * dp, struct ImplicitPolicy * pol,  */
/*                      double weight, */
/*                      int verbose, const struct ApproxArgs * aargs) */
/* { */
/*     struct DPPOL dppol; */
/*     dppol.dp = dp; */
/*     dppol.pol = pol; */
/*     dppol.weight = weight; */
/*     struct Cost * oc = dp->cost; */

/*     /\* size_t d = cost_get_d(oc); *\/ */
/*     /\* double * lb = cost_get_lb(oc); *\/ */
/*     /\* double * ub = cost_get_ub(oc); *\/ */
/*     /\* struct Cost * cost = cost_alloc(d,lb,ub); *\/ */
/*     /\* cost_init_discrete(cost,oc->N,oc->x); *\/ */

/*     struct Cost * cost = cost_copy_deep(oc); */
/* //    printf("got first one!\n"); */
/*     cost_approx(cost,dpih_rhs_pol_cost_weight,&dppol,verbose,aargs); */
/*     return cost; */
/* } */
