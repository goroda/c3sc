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

int drift_eval(struct Drift * b, double time, double * x, double * u,
               double * out)
{
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

int dyn_eval(struct Dyn * dyn, double time, double * x,
             double * u, double * drift, double * diff)
{
    int res;
    res = drift_eval(dyn->drift,time,x,u,drift);
    if (res != 0){
        return res;
    }
    res = diff_eval(dyn->diff,time,x,u,diff);
    return res;
}

size_t dyn_getdw(struct Dyn * dyn)
{
    assert (dyn != NULL);
    return diff_getdw(dyn->diff);
}

int std_dyn_eval(struct StdDyn * sd, double time,
                 double * x, double * u,
                 double * drift, double * diff)
{
    
    for (size_t ii = 0; ii < sd->d; ii++){
        sd->space[ii] = sd->slope[ii]*x[ii] + sd->off[ii];
    }
    
    int res = dyn_eval(sd->dyn,time,sd->space,u,drift,diff);
    if (res == 0){
        return res;
    }
    else{
        for (size_t ii = 0; ii < sd->d; ii++){
            drift[ii]/=sd->slope[ii];
            diff[ii]/= sd->slope[ii]; // not sure about this scaling for stochastic term;
        }
    }
    return res;
}
