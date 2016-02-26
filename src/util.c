#include <stdlib.h>
#include <stdio.h>

#include "util.h"

int c3sc_check_bounds(size_t dx, double * lbx,
                      double * ubx, double * x)
{
    if ((lbx == NULL) || (ubx == NULL)){
        return 0;
    }
    for (size_t ii = 0; ii < dx; ii++){
        if (x[ii] < lbx[ii]){
            return - (ii+1);
        }
        else if (x[ii] > ubx[ii]){
            return ii+1;
        }
    }
    return 0;
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
