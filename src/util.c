#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "util.h"
#include "cdyn/src/simulate.h"

struct ApproxArgs
{
    double cross_tol;
    double round_tol;
    size_t kickrank;
    size_t maxrank;
};

struct ApproxArgs * approx_args_init()
{
    struct ApproxArgs * aargs = malloc(sizeof(struct ApproxArgs));
    if (aargs == NULL){
        fprintf(stderr,"Failure allocating approximation arguments\n");
        exit(1);
    }
    
    aargs->cross_tol = 1e-10;
    aargs->round_tol = 1e-10;
    aargs->kickrank = 10;
    aargs->maxrank = 40;

    return aargs;
}

void approx_args_free(struct ApproxArgs * aargs)
{
    if (aargs != NULL){
        free(aargs); aargs = NULL;
    }
}

void approx_args_set_cross_tol(struct ApproxArgs * aargs, double cross_tol)
{
    assert (aargs != NULL);
    aargs->cross_tol = cross_tol;
}

double approx_args_get_cross_tol(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->cross_tol;

}
void approx_args_set_round_tol(struct ApproxArgs * aargs, double round_tol)
{
    assert (aargs != NULL);
    aargs->round_tol = round_tol;
}

double approx_args_get_round_tol(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->round_tol;
}

void approx_args_set_kickrank(struct ApproxArgs * aargs, size_t kickrank)
{
    assert (aargs != NULL);
    aargs->kickrank = kickrank;
}

size_t approx_args_get_kickrank(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->kickrank;
}

void approx_args_set_maxrank(struct ApproxArgs * aargs, size_t maxrank)
{
    assert (aargs != NULL);
    aargs->maxrank = maxrank;
}

size_t approx_args_get_maxrank(const struct ApproxArgs * aargs)
{
    assert (aargs != NULL);
    return aargs->maxrank;
}

int c3sc_check_bounds(size_t dx, double * lbx,
                      double * ubx,const  double * x)
{
    if ((lbx == NULL) || (ubx == NULL)){
        return 0;
    }
    for (size_t ii = 0; ii < dx; ii++){
        if (x[ii] < lbx[ii]){
//            printf("x-lbx = %G",x[ii]-lbx[ii]);
            return - (ii+1);
        }
        else if (x[ii] > ubx[ii]){
            return ii+1;
        }
    }
    return 0;
}

static int compareDouble (const void * a, const void * b)
{
    if (*(double*)a < *(double*)b){
        return -1;
    }
    else if (*(double*)a > *(double*)b){
        return 1;
    }
    return 0;
}

double * c3sc_combine_and_sort(size_t Nx, double * x, size_t Ny, double * y,
                               size_t * Ntot)
{
    double * temp = malloc( (Nx+Ny)*sizeof(double));
    assert (temp != NULL);
    memmove(temp,x,Nx*sizeof(double));
    memmove(temp+Nx,y,Ny*sizeof(double));
    qsort(temp,Nx+Ny,sizeof(double),compareDouble);

    double * out = calloc(Nx+Ny,sizeof(double));
    *Ntot = 1;
    out[0] = temp[0];
    for (size_t ii = 1; ii < Nx+Ny; ii++){
        if (fabs(temp[ii]-out[*Ntot-1]) > 1e-15){
            out[*Ntot] = temp[ii];
            *Ntot = *Ntot + 1;
        }
    }
    free(temp); temp = NULL;
    return out;
}


struct c3sc_SortCouple
{
    size_t ind;
    double val;
};
    
int c3sc_compare (const void * a, const void * b)
{
    struct c3sc_SortCouple *as = (struct c3sc_SortCouple *)a;
    struct c3sc_SortCouple *bs = (struct c3sc_SortCouple *)b;
    if ( as->val < bs->val ){
      return -1;  
    } 
    else if ( as->val > bs->val ){
      return 1;  
    }
    else{
        return 0;
    }
}

// probs is overwritten
size_t c3sc_sample_discrete_rv(size_t n, double * probs,
                               double sample)
{
    struct c3sc_SortCouple sc[1000];        
    if (n > 1000){
        fprintf(stderr,"Not enough memory allocate in ");
        fprintf(stderr,"discrete_sample\n");
        exit(1);
    }

    for (size_t ii = 0; ii < n; ii++){
        sc[ii].ind = ii;
        sc[ii].val = probs[ii];
//        perm[ii] = ii+1;
    }

    // sort in ascending order
    qsort (sc,n,sizeof(struct c3sc_SortCouple),c3sc_compare);

    probs[0] = sc[0].val;
    for (size_t ii = 1; ii < n; ii++){
        probs[ii] = sc[ii].val + probs[ii-1];
    }
    for (size_t jj = 0; jj < n; jj++){
        if (sample <= probs[jj]){
            return sc[jj].ind;
        }
    }
    
    fprintf(stderr,"Warning: problem with gen. sample\n");
    fprintf(stderr,"Uniform sample is %G\n",sample);
    return 0;
}
