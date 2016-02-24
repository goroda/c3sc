#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "c3.h"

#include "integrate.h"
#include "dynamics.h"
#include "simulate.h"

struct State * 
euler_step(struct State * s, struct Control * u, double dt,
           struct Dyn * dyn, double * driftv)
{
    size_t d = state_getd(s);
    double time = state_gett(s);
    double * x = state_getx_ref(s);
    double * uu = control_getu_ref(u);
    int res = dyn_eval(dyn,time,x,uu,driftv,NULL);
    if (res != 0){
        return NULL;
    }
       
    struct State * new = state_alloc();
    state_init_zero(new,d,time+dt);
    double * newx = state_getx_ref(new);
    for (size_t ii = 0; ii < d; ii++){
        newx[ii] = x[ii] + dt * driftv[ii];
    }
    
    return new;
}

struct State * 
euler_maruyama_step(struct State * s, double * noise,
                   struct Control * u, double dt,
                   struct Dyn * dyn,
                   double * drift, double * diff)
{

    size_t d = state_getd(s);
    size_t dw = dyn_getdw(dyn);
    double time = state_gett(s);
    double * x = state_getx_ref(s);
    double * uu = control_getu_ref(u);
    int res = dyn_eval(dyn,time,x,uu,drift,diff);
    if (res != 0){
        return NULL;
    }

    struct State * new = state_alloc();
    state_init_zero(new,d,time+dt);
    double * newx = state_getx_ref(new);
    for (size_t ii = 0; ii < d; ii++ ){
        double val = cblas_ddot(dw,diff+ii,d,noise,1);
        //printf("val = %G\n",val);
        newx[ii] = x[ii] + dt * drift[ii] + val;
            
    }
    return new;
}

/* int simulator_step(struct Simulator * sim, int type) */
/* { */

/*     if (type == 0){ */
/*         sim-> */
        
/*     } */
/* } */
