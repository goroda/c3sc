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

struct DPih;
{
    struct MCA * mm;
    struct Cost * cost;
    struct Policy * pol;

    double beta; // discount factor
    int (*stagecost)(double *,double *,double *);
    int (*boundcost)(double *,double *);
};

///////////////////////////////////////////////////////////////
struct DPih * dpih_alloc(double beta,
                         int (*s)(double*,double*,double*),
                         int (*b)(double*,double*))
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
}

void dpih_free(struct DPih * dp)
{
    if (dp != NULL){
        free(dp); dp = NULL;
    }
}

void dpih_attach_mca(struct DPih * dp, struct MCA * mm)
{
    assert(dp!= NULL);
    dp->mm = mca;
}

void dpih_attach_cost(struct DPih * dp, struct Cost * cost)
{
    assert(dp != NULL);
    dp->cost = cost;
}

void dpih_attach_policy(struct DPih * dp, struct Policy * pol)
{
    assert (dp != NULL);
    dp->pol = pol;
}

double dpih_rhs(struct DPih * dp,double * xin,double * u)
{

    double dt;

        int res = dp->boundcost(t,x,&out);
        assert (res == 0);
        //printf("check boundary for point ");
        //dprint(dp->bound->d, x);
        //printf("in here out=%G\n",out);
    }
    else{
        // possibly recoverable
        //int res = tensor_mm_dyn_eval(dp->tens,t,x,u);
        //double drift * tensor_mm_drift_ref(dp->tens);
        //assert (res == 0);

        // check direction of movement
        int dirin = 0;
        if (dirin == 1)
        { //direction of movement is into the domain so ok
            fprintf(stderr, "movement away from boundary not yet implemented\n")  ;
        }
        else{
            // unrecoverable
            int res = dp->boundcost(t,x,&out);
            assert (res == 0);
        }
        //printf("out = %G\n",out);
    }

    return out;
}


/* double dpih_pi_rhs_approx(double * x, void * arg) */
/* { */
/*     struct DPih * dp = arg; */
/*     struct Control * u = NULL; */
/*     double t = 0.0; */
/*     int res = policy_eval(dp->pol,t,x,&u); */
/*     assert (res == 0); */

/*     double * uu = control_getu_ref(u); */
/*     double val = dpih_rhs(dp,t,x,uu); */
    
/*     control_free(u); u = NULL; */

/*     return val; */
/* } */


/* double dpih_pi_iter_approx(struct DPih * dp,int verbose) */
/* { */

/*     struct FunctionTrain * cost =  */
/*         cost_approx_new(dp->cost,dpih_pi_rhs_approx, */
/*                         dp,verbose); */

/*     double diff = function_train_norm2diff(cost, */
/*                                            dp->cost->cost); */
/*     cost_update_ref(dp->cost,cost); */

/*     return diff; */
/* } */
