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
 *  stage cost
 *  \var DPih::boundcost
 *  boundary cost
 */
struct DPih
{
    struct MCA * mm;
    struct Cost * cost;
    struct Policy * pol;

    double beta; // discount factor
    int (*stagecost)(double,double *,double *,double *);
    int (*boundcost)(double,double *,double *);
};

///////////////////////////////////////////////////////////////

/**********************************************************//**
    Allocate an infinite horizon Dynamic program

    \param[in] beta - discount factor
    \param[in] s    - stagecost
    \param[in] b    - boundary cost

    \return dynamic program
**************************************************************/
struct DPih * dpih_alloc(double beta,
                         int (*s)(double,double*,double*,double*),
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
double dpih_rhs(struct DPih * dp,double * x,double * u)
{

    double dt;
    double t = 0.0;
    int absorb = 0;
    double val = mca_expectation(dp->mm,t,x,u,&dt,
                                 cost_eval_bb,dp->cost,
                                 &absorb);

    double out;
    if (absorb == 1){
        //printf("state absorbed\n");
        //dprint(2,x);
        int res = dp->boundcost(t,x,&out);
        //printf("cost = %G\n",out);
        assert (res == 0);
    }
    else{
        //printf("state not absorbed\n");
        //dprint(2,x);
        double sc;
        int res = dp->stagecost(t,x,u,&sc);
        assert (res == 0);
        out = exp(-dp->beta*dt)*val + dt*sc;
        //printf("dt = %G, out = %G\n",dt,out);
    }

    return out;
}

double dpih_rhs_pol(struct DPih * dp, double * x)
{
    struct Control * u = NULL;
    double t = 0.0;
    int res = policy_eval(dp->pol,t,x,&u);
    assert (res == 0);

    double * uu = control_getu_ref(u);
    double val = dpih_rhs(dp,x,uu);
    
    control_free(u); u = NULL;
    return val;    
}

double dpih_rhs_bb(double * x, void * dp)
{
    double ret = dpih_rhs_pol(dp,x);
    return ret;
}

struct Cost * dpih_iter_pol(struct DPih * dp,int verbose)
{

    struct Cost * oldcost = dp->cost;
    struct Cost * cost = cost_alloc(oldcost->d,oldcost->bds->lb,oldcost->bds->ub);
    cost_init_discrete(cost,oldcost->N,oldcost->x);
    cost_approx(cost,dpih_rhs_bb,dp,verbose);
    return cost;

}
