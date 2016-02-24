#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "dynamics.h"

void drift_init(struct Drift * b, size_t dx, size_t du,
                double * lbx, double * ubx,
                double * lbu, double * ubu)
{
    assert (b != NULL);
    b->dx = dx;
    b->du = du;
    b->lbx = lbx;
    b->ubx = ubx;
    b->lbu = lbu;
    b->ubu = ubu;
    b->b = NULL;
    b->bargs = NULL;
}

size_t drift_getdx(struct Drift * b)
{
    assert (b != NULL);
    return b->dx;
}

int drift_eval(struct Drift * b, double time, double * x, double * u,
               double * out)
{
    assert ( b != NULL);
    int res;
    if (b->b == NULL){
        fprintf(stderr,"Warning: drift dynamics (Drift->b) are not\n");
        fprintf(stderr,"         yet specified\n");
        return 1;
    }
    res = b->b(time,x,u,out,b->bargs);
    return res;
}

void diff_init(struct Diff * s, size_t dw,
               size_t dx, size_t du,
               double * lbx, double * ubx,
               double * lbu, double * ubu)
{
    assert (s != NULL);
    s->dw = dw;
    s->dx = dx;
    s->du = du;
    s->lbx = lbx;
    s->ubx = ubx;
    s->lbu = lbu;
    s->ubu = ubu;
    s->s = NULL;
    s->sargs = NULL;
}

int diff_eval(struct Diff * b, double time, double * x, 
              double * u, double * out)
{
    int res;
    if (b->s == NULL){
        fprintf(stderr,"Warning: Diff dynamics (Diff->s) are not\n");
        fprintf(stderr,"         yet specified\n");
        return 1;
    }
    res = b->s(time,x,u,out,b->sargs);
    return res;
}

size_t diff_getdw(struct Diff * b){
    assert (b != NULL);
    return b->dw;
}

void dyn_init_ref(struct Dyn *dyn, struct Drift *drift, 
                  struct Diff *diff)
{
    dyn->trans = 0;
    dyn->drift = drift;
    dyn->diff = diff;
    
    dyn->lt = NULL;
    dyn->space = NULL;
}

void dyn_add_transform_ref(struct Dyn * dyn, struct LinTransform * lt,
                           double * space)
{
    dyn->trans = 1;
    dyn->lt = lt;
    dyn->space = space;
}

size_t dyn_getdx(struct Dyn * dyn)
{
    assert (dyn != NULL);
    return drift_getdx(dyn->drift);
}

size_t dyn_getdw(struct Dyn * dyn)
{
    assert (dyn != NULL);
    return diff_getdw(dyn->diff);
}

static int dyn_eval_base(struct Dyn * dyn, double time, double * x,
                         double * u, double * drift, double * diff)
{
    int res = 0;
    if (drift != NULL){
        res = drift_eval(dyn->drift,time,x,u,drift);
    }
    if (res != 0){
        return res;
    }
    if (diff != NULL){
        res = diff_eval(dyn->diff,time,x,u,diff);
    }
    return res;
}

static int dyn_eval_std(struct Dyn * sd, double time,
                        double * x, double * u,
                        double * drift, double * diff)
{
    
    size_t dx = dyn_getdx(sd);
    lin_transform_eval(sd->lt,x,sd->space);
    int res = dyn_eval_base(sd,time,sd->space,u,drift,diff);
    if (res == 0){
        return res;
    }
    else{
        if (drift != NULL){
            for (size_t ii = 0; ii < dx; ii++){
                drift[ii] /= lin_transform_get_slopei(sd->lt,ii);
            }
            // not sure about this scaling for stochastic term;
            //diff[ii]  /= lin_transform_get_slopei(sd->lt,ii);
        }
    }
    return res;
}

int dyn_eval(struct Dyn * dyn, double time, double * x,
             double * u, double * drift, double * diff)
{
    int res;
    //printf("dyn->trans = %d\n",dyn->trans);
    if (dyn->trans == 0){
        res = dyn_eval_base(dyn,time,x,u,drift,diff);
    }
    else{
        res = dyn_eval_std(dyn,time,x,u,drift,diff);
    }
    return res;

}

