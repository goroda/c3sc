#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "c3.h"

#include "dynamics.h"

/** \struct Drift
 *  \brief Drift dynamics
 *  \var Drift::dx
 *  dimension of state space
 *  \var Drift::du
 *  dimension of control
 *  \var Drift::lbx
 *  Lower bounds for state space
 *  \var Drift::ubx
 *  Upper bounds for state space
 *  \var Drift::lbu
 *  Lower bounds for control space
 *  \var Drift::ubx
 *  Upper bounds for control space
 *  \var Drift::b
 *  RHS of drift term to stochastic differential equation
 *  f(time,state,control,out,grad,args)
 *  \var Drift::bargs
 *  Additional arguments to dynamics
 */
struct Drift
{
    size_t dx;
    size_t du;
    //these could be null
    double * lbx;
    double * ubx;
    double * lbu;
    double * ubu;

    int (*b)(double,double *, double *, double *,double*,void *);
    void * bargs;
};

struct Drift * drift_alloc(size_t dx, size_t du)
{
    struct Drift * b = malloc(sizeof(struct Drift));
    if (b == NULL){
        fprintf(stderr,"Error allocating memory for Drift\n");
        exit(1);
    }

    b->dx = dx;
    b->du = du;
    b->lbx = NULL;
    b->ubx = NULL;
    b->lbu = NULL;
    b->ubu = NULL;
    b->b = NULL;
    b->bargs = NULL;

    return b;
}

void drift_free(struct Drift * drift)
{
    if (drift != NULL){
        free(drift->lbx); drift->lbx = NULL;
        free(drift->ubx); drift->ubx = NULL;
        free(drift->lbu); drift->lbu = NULL;
        free(drift->ubu); drift->ubu = NULL;
        free(drift); drift = NULL;
    }
}

void drift_add_func(struct Drift * dr,
                    int (*b)(double,double*,double*,
                             double*,double*,void*),
                    void * bargs)
{
    assert (dr != NULL);
    dr->b = b;
    dr->bargs = bargs;
}

size_t drift_get_dx(struct Drift * b)
{
    assert (b != NULL);
    return b->dx;
}

size_t drift_get_du(struct Drift * b)
{
    assert (b != NULL);
    return b->du;
}

int drift_eval(struct Drift * b,double time,double * x,
               double * u, double * out, double * jac)
{
    assert ( b != NULL);
    int res;
    if (b->b == NULL){
        fprintf(stderr,"Warning: drift dynamics (Drift->b) are not\n");
        fprintf(stderr,"         yet specified\n");
        return 1;
    }
    res = b->b(time,x,u,out,jac,b->bargs);
    return res;
}

////////////////////////////////////////////////////////////////////////

/** \struct Diff
 *  \brief Diffusion dynamics
 *  \var Diff::dx
 *  dimension of state space
 *  \var Diff::du
 *  dimension of control
 *  \var Diff::dw
 *  dimension of random walk
 *  \var Diff::lbx
 *  Lower bounds for state space
 *  \var Diff::ubx
 *  Upper bounds for state space
 *  \var Diff::lbu
 *  Lower bounds for control space
 *  \var Diff::ubx
 *  Upper bounds for control space
 *  \var Diff::s
 *  RHS of diffusion term of stochastic differential equation
 *  f(time,state,control,out,grad,args)
 *  \var Diff::sargs
 *  Additional arguments to dynamics
 */
struct Diff
{

    size_t dx;
    size_t du;
    size_t dw;
    double * lbx;
    double * ubx;
    double * lbu;
    double * ubu;

    int (*s)(double,double*,double*,double*,double*,void*);
    void * sargs;
};


struct Diff * diff_alloc(size_t dx, size_t du, size_t dw)
{
    struct Diff * s = malloc(sizeof(struct Diff));
    if (s == NULL){
        fprintf(stderr,"Error allocating memory for Drift\n");
        exit(1);
    }

    s->dx = dx;
    s->du = du;
    s->dw = dw;
    s->lbx = NULL;
    s->ubx = NULL;
    s->lbu = NULL;
    s->ubu = NULL;
    s->s = NULL;
    s->sargs = NULL;

    return s;
}

void diff_free(struct Diff * diff)
{
    if (diff != NULL){
        free(diff->lbx); diff->lbx = NULL;
        free(diff->ubx); diff->ubx = NULL;
        free(diff->lbu); diff->lbu = NULL;
        free(diff->ubu); diff->ubu = NULL;
        free(diff); diff = NULL;
    }
}

void diff_add_func(struct Diff * df,
                   int (*s)(double,double*,double*,
                            double*,double*,void*),
                    void * sargs)
{
    assert (df != NULL);
    df->s = s;
    df->sargs = sargs;
}

int diff_eval(struct Diff * b, double time, double * x, 
              double * u, double * out, double * jac)
{
    int res;
    if (b->s == NULL){
        fprintf(stderr,"Warning: Diff dynamics (Diff->s) are not\n");
        fprintf(stderr,"         yet specified\n");
        return 1;
    }
    res = b->s(time,x,u,out,jac,b->sargs);
    return res;
}

size_t diff_get_dw(struct Diff * b){
    assert (b != NULL);
    return b->dw;
}

/////////////////////////////////////////////////////////
struct Dyn * dyn_alloc(struct Drift * drift, struct Diff * diff)
{
    struct Dyn * dyn = malloc(sizeof(struct Dyn));
    if (dyn == NULL){
        fprintf(stderr,"Cannot allocate memory for Dynamics\n");
        exit(1);
    }
    dyn->drift = drift;
    dyn->diff = diff;

    return dyn;
}

void dyn_free(struct Dyn * dyn)
{
    if (dyn != NULL){
        free(dyn); dyn = NULL;
    }
}

void dyn_free_deep(struct Dyn * dyn)
{
    if (dyn != NULL){
        drift_free(dyn->drift); dyn->drift = NULL;
        diff_free(dyn->diff); dyn->diff = NULL;
        dyn_free(dyn); dyn = NULL;
    }
}

void dyn_init_ref(struct Dyn *dyn, struct Drift *drift, 
                  struct Diff *diff)
{
    dyn->drift = drift;
    dyn->diff = diff;
}

size_t dyn_get_dx(struct Dyn * dyn)
{
    assert (dyn != NULL);
    return drift_get_dx(dyn->drift);
}

size_t dyn_get_dw(struct Dyn * dyn)
{
    assert (dyn != NULL);
    return diff_get_dw(dyn->diff);
}

size_t dyn_get_du(struct Dyn * dyn)
{
    assert (dyn != NULL);
    return drift_get_du(dyn->drift);
}

int dyn_eval(struct Dyn * dyn,double time,
             double * x,
             double * u, double * drift, double * jacdr,
             double * diff, double * jacdiff)
{
    int res = 0;
    if (drift != NULL){
        res = drift_eval(dyn->drift,time,x,u,drift,jacdr);
    }
    if (res != 0){
        return res;
    }
    if (diff != NULL){
        res = diff_eval(dyn->diff,time,x,u,diff,jacdiff);
    }
    return res;
}


