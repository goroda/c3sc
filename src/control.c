#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "c3.h"
#include "control.h"
#include "simulate.h"
#include "util.h"

struct Policy * policy_alloc()
{
    struct Policy * pol = NULL;
    pol = malloc(sizeof(struct Policy));
    if (pol == NULL){
        fprintf(stderr, "Cannot allocate Policy \n");
        exit(1);
    }
    pol->dx = 0;
    pol->du = 0.0;
    pol->lbx = NULL;
    pol->ubx = NULL;
    pol->feedback = NULL;
    pol->trans = 0;
    pol->lt = NULL;
    pol->space = NULL;
    
    return pol;
}

void policy_init(struct Policy * pol, size_t dx, size_t du, 
                 double * lbx, double * ubx)
{
    pol->dx = dx;
    pol->du = du;
    if (lbx != NULL){
        pol->lbx = calloc_double(dx);
        memmove(pol->lbx,lbx,dx*sizeof(double));
    }
    if (ubx != NULL){
        pol->ubx = calloc_double(dx);
        memmove(pol->ubx,ubx,dx*sizeof(double));
    }
   
}

void policy_add_transform_ref(struct Policy * pol,
                              struct LinTransform * lt, 
                              double * space)
{
    pol->trans = 1;
    pol->lt = lt;
    pol->space = space;
}

void policy_free(struct Policy * pol)
{
    if (pol != NULL){
        free(pol->lbx); pol->lbx = NULL;
        free(pol->ubx); pol->ubx = NULL;
        free(pol); pol = NULL;
    }
}

void policy_add_feedback(struct Policy * pol,
                          int (*f)(double,double*,double*))
{
    assert (pol != NULL);
    pol->feedback = f;
}

int policy_eval(struct Policy * pol, double time,
                double * xin, 
                struct Control ** u)
{
    assert (pol != NULL);
    assert (pol->feedback != NULL);
    double * x;
    if (pol->trans == 0){
        x = xin;
    }
    else{
        lin_transform_eval(pol->lt,xin,pol->space);
        x = pol->space;
    }
    
    int res = c3sc_check_bounds(pol->dx,pol->lbx,pol->ubx,x);
    if (res != 0){
        return res;
    }
    
    assert (*u == NULL);
    (*u) = control_alloc();
    (*u)->d = pol->du;
    (*u)->u = calloc_double(pol->du);
    res = pol->feedback(time,x,(*u)->u);
    return res;
}
