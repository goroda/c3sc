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
    dp_param_add_boundary(c3c->bound);

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
    do_param_add_drift(c3c->dp,b,args);
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
    dp_param_add_boundcost(c3c->dp,obscost);
}

