#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "c3.h"

int euler_maruyama_step(double * x, double * noise, double dt,
                        double t, double * x, double * u,
                        double * drift,double * diff,
                        struct Dyn * f)
{

    int res; 
    res = dyn_eval(f,t,x,u,drift,diff);
    if (res != 0){
        return res;
    }
    size_t dim = dyn_get_dx(f);
    size_t dimw = dyn_get_dw(f);
    
    for (size_t ii = 0; ii < dim; ii++ ){
        x[ii] = x[ii] + dt * drift[ii] + cblas_ddot(dimw,diff+ii,dim,noise,1);
    }
    
}

int simulator_step(struct Simulator * sim, int type)
{

    if (type == 0){
        sim->
        
    }
}
