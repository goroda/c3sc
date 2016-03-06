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
 */
struct DPih
{
    struct MCA * mm;
    struct Cost * cost;
    struct Policy * pol;

    double beta; // discount factor
    int (*stagecost)(double,double*,double*,double*,double*);
    int (*boundcost)(double,double*,double*);

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
        du = mca_du(dp->mm);
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
            du = mca_du(dp->mm);
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
            free(gdt); gdt = NULL;
        }

    }

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
    
    size_t du = mca_du(dpx.dp->mm);
    double * ustart = calloc_double(du);
    // need to include boundaries
    /* double * lb = cost_get_lb(dpx.dp->cost); */
    /* double * ub = cost_get_ub(dpx.dp->cost); */
    
    double val = 0.0;
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    /* c3opt_add_lb(opt,lb); */
    /* c3opt_add_ub(opt,ub); */
    c3opt_set_verbose(opt,0);
    c3opt_add_objective(opt,dpih_rhs_opt_bb,&dpx);
    


    int res = c3opt_minimize(opt,ustart,&val);
    /* if (res != 0){ */
    /*     dprint(2,x); */
    /*     dprint(1,ustart); */
    /*     printf("cost = %G\n",val); */
    /* } */
    /* assert (res == 0); */
    //double val = OPTIMIZE(du,ustart,dpih_rhs_opt_bb,&dpx);
    
    free(ustart); ustart = NULL;
    c3opt_free(opt); opt = NULL;
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
