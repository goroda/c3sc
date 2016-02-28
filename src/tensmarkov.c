#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "c3.h"

#include "tensmarkov.h"
#include "util.h"
#include "simulate.h"
#include "dynamics.h"
#include "cost.h"

void tensor_mm_init_ref(struct TensorMM *tens, size_t d, 
                        double h, struct Dyn * dyn, double * space)
{
    tens->d = d;
    tens->h = h;
    tens->dyn = dyn;
    tens->space = space;
}

int tensor_mm_tprob(struct TensorMM * tens, double t,
                    double * x, double * u,
                    double * probs, double * dt)
{
    // NOTE ONLY WORKS for diagonal based diffusions
    // probs are +- state for each dimension
    size_t dim = tens->d;
    double h = tens->h;
    //printf("eval dyn\n");
    int res = dyn_eval(tens->dyn,t,x,u,tens->space,
                       tens->space+dim);
    //printf("res = %d\n",res);
    if (res != 0){
        return res;
    }
    
    size_t dw = dyn_getdw(tens->dyn);
    double norm = 0;
    double diff;
    for (size_t ii = 0; ii < dim; ii++){
        diff = cblas_ddot(dw,
                          tens->space+dim+ii,dim,
                          tens->space+dim+ii,dim);
        probs[ii*2] = diff/2.0;
        probs[ii*2+1] = diff/2.0;
        norm += diff;
        if (tens->space[ii] > 0){
            probs[ii*2+1] += h * tens->space[ii];
            norm += h * tens->space[ii];
        }
        else{
            probs[ii*2] += h *(-tens->space[ii]);
            norm += h * (-tens->space[ii]);
        }
    }
    *dt = h*h / norm;
    probs[dim*2] = 1.0;
    for (size_t ii = 0; ii < dim ; ii++){
        probs[ii*2] /= norm;
        probs[ii*2+1] /= norm;
        probs[dim*2] -= probs[ii*2];
        probs[dim*2] -= probs[ii*2+1];
    }
    return 0;
}

struct State *
tensor_mm_step(struct TensorMM * mm,
               struct State * s, double noise,
               struct Control * u, double * dt,
               double * probs)
{
    size_t d = dyn_getdx(mm->dyn);
    double time = state_gett(s);
    double * x = state_getx_ref(s);
    double * uu = control_getu_ref(u);
//    printf("state = "); dprint(d,x);
    int res = tensor_mm_tprob(mm,time,x,uu,probs,dt);
//    printf("probs are = "); dprint(2*d+1,probs);
//    printf("got probs\n");

    if (res != 0){
        return NULL;
    }

    size_t ind = c3sc_sample_discrete_rv(2*d+1,probs,noise);

    struct State * new = state_alloc();
    state_init_zero(new,d,time + *dt);
    double * newx = state_getx_ref(new);
    memmove(newx,x,d*sizeof(double));
    if (ind != 2*d){
        if ((ind % 2) == 0){ // even move backwards
            newx[ind/2] = x[ind/2] - mm->h;
        }
        else{
            newx[(ind+1)/2 -1] = x[(ind+1)/2-1] + mm->h;
        }
    }
    
    return new;
    
}


// assumes only transitions to neighbors
double tensor_mm_cost(struct TensorMM * mm,
                      double t, double * x,
                      double * uu,
                      struct Cost * cost,
                      double * probs)
{
    double dt;
    int res = tensor_mm_tprob(mm,t,x,uu,probs,&dt);
    assert (res != 0);
    double evals[2];
    double pt[2];
    double out = 0.0;

    //evaluate transition to oneself
    res = cost_eval(cost,t,x,&out);
    assert (res != 0);
    out *= probs[2*mm->d];

    for (size_t ii = 0; ii < mm->h; ii++){
        pt[0] = x[ii]-mm->h;
        pt[1] = x[ii]+mm->h;
        res = cost_eval_neigh(cost,t,x,ii,pt,evals);
        assert (res != 0);
        out += probs[ii*mm->d]*evals[0];
        out += probs[ii*mm->d+1]*evals[1];
    }


    return out;
}
