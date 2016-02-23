#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "integrate.h"
#include "dynamics.h"
#include "simulate.h"

struct State * 
euler_step(struct State * s, struct Control * u, double dt,
           struct Drift * drift, double * driftv)
{
    size_t d = state_getd(s);
    double time = state_gett(s);
    double * x = state_getx_ref(s);
    double * uu = control_getu_ref(u);
    int res = drift_eval(drift,time,x,uu,driftv);
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


/* int euler_maruyama_step(double * x, double * noise, double dt, */
/*                         double t, double * u, */
/*                         double * drift,double * diff, */
/*                         struct Dyn * f) */
/* { */

/*     int res; */
/*     res = dyn_eval(f,t,x,u,drift,diff); */
/*     if (res != 0){ */
/*         return res; */
/*     } */
/*     size_t dim = dyn_get_dx(f); */
/*     size_t dimw = dyn_get_dw(f); */
    
/*     for (size_t ii = 0; ii < dim; ii++ ){ */
/*         x[ii] = x[ii] + dt * drift[ii] + cblas_ddot(dimw,diff+ii,dim,noise,1); */
/*     } */
    
/* } */

/* int simulator_step(struct Simulator * sim, int type) */
/* { */

/*     if (type == 0){ */
/*         sim-> */
        
/*     } */
/* } */
